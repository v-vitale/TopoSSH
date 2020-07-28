clear
close all
load data64.mat

%cm=[linspace(1,1,50)  linspace(1,.1,50); linspace(.1,1,50) linspace(1,.1,50); linspace(.1,1,50) linspace(1,1,50)]';
threshold=0.1;

tv=[1e-3 0.1:0.1:t(end)];
Wv= 0:0.2:W(end);

Sm=zeros(length(Wv),length(tv));
dSm=Sm;
rel=Sm;
for i=1:length(S)
    Sm(i)=S(i);
    dSm(i)=dS(i);
    if Sm(i)>threshold
        rel(i)=dSm(i)/Sm(i);
    else
        rel(i)=0;
    end
end

[Y,X]=meshgrid(tv,Wv);
f=figure;
surf(X,Y,Sm)
yticks([0 .5 1 1.5])

cmap=[linspace(1,1,500) linspace(1,40/255,500); linspace(1,102/255,500) linspace(102/255,70/255,500); linspace(0,0,500) linspace(0,235/255,500)]';
cmp=colormap(cmap);
grid off
cmp=flipud(cmp);
colormap(cmp);
h=colorbar;
xlim([0,5.8])
ylim([0,1.9])
view(2);

ax=gca;
set(gca,'linewidth',2)
box on
ax.FontSize = 20; 
xlabel('$W/w$','Interpreter','latex','Fontsize',30);
ylabel('$v/w$','Interpreter','latex','Fontsize',30);

ylabel(h,'$S^\mathrm{D}/\log 2$','Interpreter','latex','Fontsize',30);
set(h,'YTick',[0 .5 1 1.5 2])


f.Renderer='Painters';
set(f,'Units','Inches');
pos = get(f,'Position');
set(f,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
%print(f,'phase_diagram','-dpdf','-r0')
%print(gcf, 'phase_diagram.pdf', '-dpdf', '-loose')


% figure
% surf(X,Y,dSm)
% %shading interp
% colormap('jet')
% colorbar
% view(2);  
% title('standard deviation')
% 
% figure
% surf(X,Y,rel)
% %shading interp
% colormap('jet')
% colorbar
% view(2);  
% title('relative error')