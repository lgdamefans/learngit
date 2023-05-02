op1 = [92.9816687698600;94.6097342157858;96.0907268881011;97.4926343911639;98.9482877413952;100.535060459018;102.162062263852;103.724981903480;105.282227671363;106.850141424814;108.339352182345;109.843467202030;111.357761543621];
op1 = op1.';
op2 = [114.415348681846;114.415348681846;114.410920936788;114.407755024705;114.407755024705;114.402057917804;114.399361270598;114.396195358514;114.395206896637;114.395206896637;114.389509789737;114.386813142530;114.386813142530];
op2 = op2.';
pmin = 0.3:0.05:0.9;
p = plot(pmin,op1,'--',pmin,op2,'-','marker','o');
grid on;
set(gca,'fontsize',14);
xlabel('Data accuracy threshold','fontsize',14);
ylabel('System cost','fontsize',14);
ylim([85,130]);
xlim([0.25,0.95]);
legend('Proposed','Alldata','location','northwest','fontsize',14);
