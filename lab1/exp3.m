%parameters
timeslot = 10;
device = 20;
t = 1;
obj = zeros(14,4);
distance = 10+rand(1,50)*40;
geoDistance = 1+rand(1,50)*100;
randompolicy = rand(timeslot,50);
for i=1:timeslot
    for j = 1:50
        if(randompolicy(i,j)>0.5)
            randompolicy(i,j) = 1;
        else
            randompolicy(i,j) = 0;
        end
    end
end
totalGeo = zeros(timeslot,50);
totalDistance = zeros(timeslot,50);
for i = 1:timeslot
    totalDistance(i,:) = distance;
    totalGeo(i,:) = geoDistance;
end
pmin = 0.6;
for accuracy = 0.3:0.05:0.9;
totalnum = timeslot*device;
dataSize = ones(timeslot,device)*3000000; %kbits
cpuCycle = ones(timeslot,device)*1e09;%Megacycles
sensingFactor = 0.06;



sensingPos = zeros(timeslot,device);
sensingPro = zeros(timeslot,device);
sensingtime = ones(timeslot,device);
senseUnit = 0.2;
sensinglatency = sensingtime*0.2;
bandwidth = 50e06;
baseNoise = 9.999999999999962e-14;
tranPower = 0.2;
gain = zeros(timeslot,device);
translatency = zeros(timeslot,device);
edgelatency = zeros(timeslot,device);
senseEnergy = zeros(timeslot,device);
energyUnit = 1e-9;
Tmax = 3;
Emax = 0.04;
weightT = 0.5;
weightE  = 0.5;
for i = 1:timeslot
    for j = 1:device
    sensingPos(i,j) = exp(-sensingFactor*totalDistance(i,j));
    sensingPro(i,j) = 1-(1-sensingPos(i,j))^sensingtime(i,j);
    gain(i,j) = 128.1+37.6*log10(totalGeo(i,j)/1000);
    gain(i,j) = 1/(10^(gain(i,j)/10));
    translatency(i,j) = dataSize(i,j)/bandwidth*log2(1+(gain(i,j)*tranPower/baseNoise));
    edgelatency(i,j)= cpuCycle(i,j)/2e09;
    energy_edge = edgelatency.*tranPower;
    end
end
data = dataSize;

leastsensingtime = ones(timeslot,device);
for i = 1:timeslot
    for k = 1:device
    for j = 1:100
          if((1-(1-sensingPos(i,k))^j)>pmin)
            leastsensingtime(i,k) = j;
            break;
        end
      
    end
    end
end
policy = ones(timeslot,device)*0.5;
for iter = 1:3
  for i = 1:timeslot
    for j = 1:device
    sensingPro(i,j) = 1-(1-sensingPos(i,j))^sensingtime(i,j);
    translatency(i,j) = data(i,j)/bandwidth*log2(1+(gain(i,j)*tranPower/baseNoise));
    end
  end  
for i = 1:timeslot
    for k = 1:device
    pre = 10000;
    for j = 1:30
        p1 = (j*senseUnit+edgelatency(i,k)+translatency(i,k));
        p2 = ((1-(1-sensingPos(i,k))^j)*Tmax);
        p3 = (j*energyUnit*dataSize(i,k)+energy_edge(i,k))/((1-(1-sensingPos(i,k))^j)*Emax)*weightE;
        cur = (j*senseUnit+edgelatency(i,k)+translatency(i,k))/((1-(1-sensingPos(i,k))^j)*Tmax)*weightT+(j*energyUnit*dataSize(i,k)+energy_edge(i,k))/((1-(1-sensingPos(i,k))^j)*Emax)*weightE;
       disp(cur);
       disp(pre);
        if(cur>pre)
               sensingtime(i,k) = max(leastsensingtime(i,k),j-1);
               break;
        end
        
        pre = cur;
    end
    end
end
% 
data = data/10000;

datasmallsize = ones(timeslot,device)*300;            
cvx_begin
      variable data(timeslot,device)
      minimize sum(sum(data))
      subject to 
        for i = 1:timeslot
            for j = 1:device
                data(i,j)<=datasmallsize(i,j);
                data(i,j)>=0;
            end
        end
        for i = 1:timeslot
            transdata = 0;
            alldata = 0;
            for j = 1:i
                for k = 1:device
                    transdata = transdata + data(j,k)*policy(j,k);
                    alldata = alldata + datasmallsize(j,k)*policy(j,k);
                end
            end
            if(alldata>0)
                (transdata/alldata) >= accuracy;
            end
        end
 cvx_end     
 cost = zeros(timeslot,device);
 data = data*10000;
 for i = 1:timeslot
     for j = 1:device
         cost(i,j) = ((sensingtime(i,j)*senseUnit+(data(i,j)/bandwidth*log2(1+(gain(i,j)*tranPower/baseNoise)))+edgelatency(i,j))/(((1-(1-sensingPos(i,j))^sensingtime(i,j))*Tmax)))*weightT+((sensingtime(i,j)*energyUnit*data(i,j)+energy_edge(i,j))/((1-(1-sensingPos(i,j))^sensingtime(i,j))*Emax))*weightE;
     end
 end
 block = 45;
 limit = 1;
 cvx_begin
    variable x(timeslot,device)
    minimize sum(sum(cost.*x))
    for i = 1:timeslot
        for j = 1:device
            x(i,j)>=0;
            x(i,j)<=1;
        end
    end
    for i = 1:timeslot
        sum(x(i,:))<=block;
        sum(x(i,:))>0;
    end
    for j = 1:timeslot-3
        for i = 1:device
            sum(x(j:j+3,i))>=limit;
        end
    end
    for i = 1:timeslot
            transdata = 0;
            alldata = 0;
            for j = 1:i
                for k = 1:device
                    transdata = transdata + data(j,k)/10000*x(j,k);
                    alldata = alldata + dataSize(j,k)/10000*x(j,k);
                end
            end
           transdata >= accuracy*alldata;
          
    end
cvx_end
for i = 1:timeslot
        for j = 1:device
            if(x(i,j)>0.5)
                policy(i,j) = 1;
            else
                policy(i,j) = 0;
            end
        end
end
end
allcost = 0;
alldata = 0;
allsense = 0;
allecost = 0;
for i = 1:timeslot
     for j = 1:device
         cost(i,j) = ((sensingtime(i,j)*senseUnit+(data(i,j)/bandwidth*log2(1+(gain(i,j)*tranPower/baseNoise)))+edgelatency(i,j))/(((1-(1-sensingPos(i,j))^sensingtime(i,j))*Tmax)))*weightT+((sensingtime(i,j)*energyUnit*data(i,j)+energy_edge(i,j))/((1-(1-sensingPos(i,j))^sensingtime(i,j))*Emax))*weightE;
         datacost = ((sensingtime(i,j)*senseUnit+(dataSize(i,j)/bandwidth*log2(1+(gain(i,j)*tranPower/baseNoise)))+edgelatency(i,j))/(((1-(1-sensingPos(i,j))^sensingtime(i,j))*Tmax)))*weightT+((sensingtime(i,j)*energyUnit*dataSize(i,j)+energy_edge(i,j))/((1-(1-sensingPos(i,j))^sensingtime(i,j))*Emax))*weightE;
         sensecost = ((leastsensingtime(i,j)*senseUnit+(data(i,j)/bandwidth*log2(1+(gain(i,j)*tranPower/baseNoise)))+edgelatency(i,j))/(((1-(1-sensingPos(i,j))^leastsensingtime(i,j))*Tmax)))*weightT+((leastsensingtime(i,j)*energyUnit*data(i,j)+energy_edge(i,j))/((1-(1-sensingPos(i,j))^leastsensingtime(i,j))*Emax))*weightE;
         ecost = ((sensingtime(i,j)*senseUnit+(data(i,j)/bandwidth*log2(1+(gain(i,j)*tranPower/baseNoise)))+edgelatency(i,j))/(((1-(1-sensingPos(i,j))^sensingtime(i,j))*Tmax)))*weightT+((sensingtime(i,j)*energyUnit*data(i,j)+energy_edge(i,j))/((1-(1-sensingPos(i,j))^sensingtime(i,j))*Emax))*weightE;
         allecost = allecost+cost(i,j)*policy(i,j);
         allcost = allcost+cost(i,j)*randompolicy(i,j);
         alldata = alldata+datacost*policy(i,j);
         allsense = allsense+sensecost*policy(i,j);
     end
end
optimal = sum(sum(cost.*policy));
obj(t,1) = allecost;
obj(t,2) = allcost;
obj(t,3) = alldata;
obj(t,4) = allsense;
t = t+1;
end

%end
% cvx_begin
%     variable s(1,totalnum)
%     minimize sum(((s*senseUnit+translatency+edgelatency).*(1/((1-(1-sensingPos).^s)*Tmax))*weightT+(s*energyUnit.*dataSize+energy_edge).*(1/((1-(1-sensingPos).^s)*Emax))*weightE));
%     subject to
%        (1-(1-sensingPos).^s) >= pmin
%       s*senseUnit+translatency+edgelatency<=Tmax
%         s*energyUnit.*dataSize+energy_edge<=Emax
% cvx_end


