close all

v1=.1;
v2=1.5;

Lvec=[32 64 80 96 112 128 144 160 ];%192 208 224 240 256];

f=figure;
ax=gca;
set(gca,'linewidth',2)
box on
ax.FontSize = 20; 
hold all
%set(0,'DefaultAxesColorOrder',hsv(9));


for ii=1:length(Lvec)
    L=Lvec(ii);
    place=['Sq_quench_L=' num2str(L) '_v1=' num2str(v1) '_v2=' num2str(v2) '.dat'];
    M=load(place);
    plot(M(:,1),M(:,6),'LineWidth',2.5,'DisplayName',['$L=$' num2str(L)]);
end
ylim([-8,35])
ylabel('$S^\mathrm{D}/\log2$','Interpreter','latex','Fontsize',30);
xlabel('$t$','Interpreter','latex','Fontsize',30);
yticks([0 10 20 30])
legend('interpreter','latex','location','s','FontSize',21,'NumColumns',4)
pos = get(gca, 'Position');
set(gca, 'Position', [pos(1)+.001 pos(2) pos(3) pos(4)]);

axin=axes('position',[.23 .605 .285 .270]);
%axin=gca;
set(gca,'linewidth',2)
box on
axin.FontSize = 18; 
hold all
p=[30 50 70 100 100 100 120 140 ]% 100 100 100 100 100 100 100 100 100 100 100 100 100 100];
for ii=1:length(Lvec)
    L=Lvec(ii);
    place=['Sq_quench_L=' num2str(L) '_v1=' num2str(v1) '_v2=' num2str(v2) '.dat'];
    M=load(place);
    plot(M(1:p(ii),1),M(1:p(ii),6),'LineWidth',3,'DisplayName',['$L=$' num2str(L)]);
end
plot([0 20],[2 2],'-.k','LineWidth',2.5)
plot([0 20],[2.02 2.02],'-.k','LineWidth',2.5)
ylim([1.99 2.03])
xlim([0 12])

set(f,'Units','Inches');
pos = get(f,'Position');
set(f,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
%print(f,'post_quench_evol','-dpdf','-r0')