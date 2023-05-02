op1 = [124.346740675060;115.343934268455;107.030111841822;102.647103567230;102.625134889850;102.608877748380];
op1 = op1.';
%op2 = [203.842707203157,232.344913463887,275.150205231889,334.677399977851,364.816703791573,414.982351726124,438.977822617271];
op2 = [221.150051239766;220.052311984374;220.300160181396;219.866991992529;219.629868995281;219.474485348781];
op2 = op2.';
op3 = [127.092208800826;117.662353321895;109.289708074742;105.359357394520;105.359357394520;105.359357394520];
op3 = op3.';
op4 = [143.414645071846;132.859108871635;122.805859771876;117.379869280648;117.356857838087;117.339434667383];
op4 = op4.';
total = [op1;op2;op3;op4];
total = total.';
device = 5:2:15;

p = bar(device,total,1);
set(p(1),'FaceColor',[78,98,171]/255)
set(p(2),'FaceColor',[70,158,180]/255)
set(p(3),'FaceColor',[135,207,164]/255)
set(p(4),'FaceColor',[254,232,154]/255)
grid on;
set(gca,'fontsize',14);
xlabel('Number of resource blocks','fontsize',14);
ylabel('System cost','fontsize',14);
legend('Proposed','Random','Firstfit','Greedysense','location','northwest','fontsize',14);