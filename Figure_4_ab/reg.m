clear 
clc

v1=.1;
v2=1.5;

Lvec=[144 160 176 192 208 224 240 256]';

f=figure;
ax=gca;
set(gca,'linewidth',2)
box on
ax.FontSize = 20; 
hold all
% set(0,'DefaultAxesColorOrder',hsv(14));

tc=zeros(size(Lvec));
for ii=1:length(Lvec)
    L=Lvec(ii);
    place=['Sq_quench_L=' num2str(L) '_v1=' num2str(v1) '_v2=' num2str(v2) '.dat'];
    M=load(place);
    %plot(M(:,1),M(:,6),'LineWidth',2.5,'DisplayName',['$L=$' num2str(L)]);
    tc(ii)=findpoint(M(:,1),M(:,6),2.02);
end
Linv=1./Lvec(2:end);
tinv=1./tc(2:end);

cf = fit(Linv,tinv,'poly1','Lower',[-Inf 0]); 
x=linspace(-1,83e-3,20);
cf_coeff = coeffvalues(cf);
cf_confint = confint(cf);
a = cf_coeff(1);
b = cf_coeff(2);
a_uncert = (cf_confint(2,1) - cf_confint(1,1))/2;
b_uncert = (cf_confint(2,2) - cf_confint(1,2))/2;
y=a*x+b;
h1=plot(1./Lvec,1./tc,'d','MarkerSize',12,'LineWidth',3);
plot(x,y,'LineWidth',2.5);

uistack(h1,'top');
set(h1, 'markerfacecolor', get(h1, 'color'))
ylim([0 .12])
xlim([0 8e-3])
ylabel('$1/t_\mathrm{c}$','Interpreter','latex','Fontsize',30);
xlabel('$1/L$','Interpreter','latex','Fontsize',30);
%yticks([0 10 20 30])
%legend('interpreter','latex','location','nw','FontSize',21)
pos = get(gca, 'Position');
set(gca, 'Position', [pos(1)+.001 pos(2) pos(3) pos(4)]);



set(f,'Units','Inches');
pos = get(f,'Position');
set(f,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
print(f,'tc_scaling','-dpdf','-r0')
function out=findpoint(x,y,p)
    r=y-p;
    s=sign(r);
    d=diff(s);
    ind=find(d~=0,1);
    y0=r(ind);
    y1=r(ind+1);
    x0=x(ind);
    x1=x(ind+1);
    
    out=x0-y0*(x1-x0)/(y1-y0);
end