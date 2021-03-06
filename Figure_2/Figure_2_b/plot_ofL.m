close all
misizebig=8;
msizesmall=4;
f=figure;
COL=get(gca,'colororder');
ax=gca;
box on
ax.FontSize = 16; 

M=load('SD_of_L_r=0.1.dat');
hold on
plot(M(:,1), 2-M(:,2),'s',...
    'MarkerFaceColor',COL(1,:),...
    'MarkerEdgeColor',COL(1,:),...
    'MarkerSize',misizebig,...
    'DisplayName','$v/w=0.1$')
plot(M(:,1), 2-M(:,3),'s-',...
    'Color',COL(1,:),...
    'LineWidth',2,...
    'MarkerFaceColor',COL(1,:),...
    'MarkerEdgeColor',COL(1,:),...
    'MarkerSize',msizesmall,...
    'HandleVisibility','off')

M=load('SD_of_L_r=0.5.dat');
plot(M(:,1), 2-M(:,2),'o',...
    'MarkerFaceColor',COL(2,:),...
    'MarkerEdgeColor',COL(2,:),...
    'MarkerSize',misizebig,...
    'DisplayName','$v/w=0.5$')
plot(M(:,1), 2-M(:,3),'o-',...
    'Color',COL(2,:),...
    'LineWidth',2,...
    'MarkerFaceColor',COL(2,:),...
    'MarkerEdgeColor',COL(2,:),...
    'MarkerSize',msizesmall,...
    'HandleVisibility','off')

M=load('SD_of_L_r=2.dat');
plot(M(:,1), M(:,2),'d',...
    'MarkerFaceColor',COL(3,:),...
    'MarkerEdgeColor',COL(3,:),...
    'MarkerSize',misizebig,...
    'DisplayName','$v/w=2$')

M=load('SD_of_L_r=10.dat');
plot(M(:,1), M(:,2),'v',...
    'MarkerFaceColor',COL(4,:),...
    'MarkerEdgeColor',COL(4,:),...
    'MarkerSize',misizebig,...
    'DisplayName','$v/w=10$')
set(gca,'linewidth',2)
ax.YAxis(1).Scale='log';
ax.YAxis(1).MinorTick='off';

yl=ylabel('$2-S^\mathrm{D}/\log2$','Interpreter','latex','Fontsize',30);
%yl.Position(1) = yl.Position(1) + 0.001;
xlabel('$L$','Interpreter','latex','Fontsize',30);

ylim([1e-25,10])
set(ax,'YTick',[1e-24,1e-16,1e-9,1e-2]);
set(ax,'XTick',[16,32,48,64]);
pos = get(gca, 'Position');
set(gca, 'Position', [pos(1)+.01 pos(2) pos(3) pos(4)]);
%right left axis
yyaxis right
ax.YAxis(2).Color='k';
ax.YAxis(2).Scale='log';
ax.YAxis(2).MinorTick='off';
ylabel('$S^\mathrm{D}/\log2$','Interpreter','latex','Fontsize',30);


ylim([1e-25,10])
set(ax,'YTick',[1e-24,1e-16,1e-9,1e-2]);
set(gca,'Yticklabel',[]) 


legend('interpreter','latex','location','sw','FontSize',29)

set(f,'Units','Inches');
pos = get(f,'Position');
set(f,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
print(f,'SD_of_L','-dpdf','-r0')