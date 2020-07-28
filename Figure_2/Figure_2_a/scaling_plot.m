clear
clc
close all
mymarkers='odsvp';
mml=length(mymarkers);
Lv=16:32:160;
f=figure;
ax=gca;
set(gca,'linewidth',2)
box on
ax.FontSize = 20; 
hold all
for i=1:length(Lv)
    L=Lv(i);
    place=['Sq_function_of_hopping_L=' num2str(L) '.dat'];
    M=load(place);
    h=plot(M(:,1),M(:,2),'--','DisplayName',['$L=$' num2str(L)]);
    set(h,'Marker',mymarkers(mod(i-1,mml)+1))
    set(h,'MarkerSize',7)
    set(h,'LineWidth',2)
    set(h, 'MarkerFaceColor', get(h,'Color')); 
end
set(ax,'YTick',[0 1 2]);

ylabel('$S^\mathrm{D}/\log2$','Interpreter','latex','Fontsize',30);
ylim([-.1,2.1])

xlabel('$v/w$','Interpreter','latex','Fontsize',30);

legend('interpreter','latex','location','sw','FontSize',23)

set(f,'Units','Inches');
pos = get(f,'Position');
set(f,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
print(f,'scaling_even','-dpdf','-r0')