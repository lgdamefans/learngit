%
total = [116.543306406895,229.909754964333,119.787888264111,127.812398239868;139.552292407468,273.678414314081,143.337903093228,154.842642370178;166.716897245278,344.318145587933,171.220036510720,185.346811750215;195.198185740003,415.827178644694,200.483788476112,216.956915338211;221.343054705909,469.395313883007,227.364046227249,246.088515110942;246.713652682674,530.597371675859,253.315521657681,274.627484718140;277.942704846061,581.386122517560,285.487140149094,308.099202158893];
device = [20,25,30,35,40,45,50];

p = bar(device,total,1);
set(p(1),'FaceColor',[78,98,171]/255)
set(p(2),'FaceColor',[70,158,180]/255)
set(p(3),'FaceColor',[135,207,164]/255)
set(p(4),'FaceColor',[254,232,154]/255)
grid on;
set(gca,'fontsize',14);
xlabel('Number of physical devices','fontsize',14);
ylabel('System cost','fontsize',14);
legend('Proposed','Random','Firstfit','Greedysense','location','northwest','fontsize',14);