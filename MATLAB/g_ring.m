%% Electrostatic reconstruction of secular equation for ring+g
D=1;
L=1;
G=1e-3;
N=1000;
g=G/N;
epsilon = linspace(0,5000,1e5);
s=20;
S=s*L;

q = sqrt(epsilon/D-s^2/4)*L;
epsilon_1 =0; S^2/4 * G/(G-1)
% epsilon_1=0.0635;1.5*1e-5;0.25
psi =( 2*(cos(q)+ 1/G * (q.^2+S^2/4)./(2*q) .*sin(q))-2);
V = log(2*(sqrt(1+(1/G * (q.^2+S^2/4)./(4*q)).^2)-2));


%find epsilon_k from continuuous ring and calculate electrostatic potential
ind = find(abs(diff(sign(psi)))==2); 
ind2 = ind(1:10);
% [kappa,signKappa,expKappa,E_vec] = thoulessKappa([epsilon_1 epsilon(ind)],N,g,[ epsilon]);
% [kappa2,signKappa2,expKappa2,E_vec2] = thoulessKappa([epsilon_1 epsilon(ind2)],N,g,[ epsilon]);
[kappa,signKappa,expKappa,E_vec] = thoulessKappa([0 epsilon(ind)],N,g,[ epsilon]);
[kappa2,signKappa2,expKappa2,E_vec2] = thoulessKappa([0 epsilon(ind2)],N,g,[ epsilon]);

[kappa_wrong,signKappa_wrong,expKappa_wrong,E_vec] = thoulessKappa([ epsilon(ind)],N,g,[ epsilon]);

figure;
axes('FontSize',24)
hold on
% plot(epsilon,psi,epsilon,S/2*ones(1,length(epsilon)),epsilon(ind),zeros(1,length(ind)),'*')
h1 = plot(epsilon,log(psi),'b','LineWidth',4);...
% h3 = plot(epsilon(ind),zeros(1,length(ind)),'*k','LineWidth',4)
% plot(epsilon,V)
% ax=get(gca);
% plot(E_vec,-E_vec'.*expKappa/g,'r')
h2 = plot(E_vec,((kappa)-kappa(2)+ log(2*cosh(S/2)-2)),'g','LineWidth',4);
h4 = plot(E_vec,((kappa_wrong)-kappa_wrong(2)+ log(2*cosh(S/2)-2)),'y','LineWidth',4)

h3=plot(E_vec,((kappa2)-kappa2(2)+ log(2*cosh(S/2)-2)),'--k','LineWidth',2)
plot(epsilon,log(2*cosh(S/2)-2)*ones(1,length(epsilon)),'r--','LineWidth',4)


% plot(E_vec,kappa,'r')
% plot(E_vec,(log(sqrt(expKappa))),'g')
% plot([epsilon_1,epsilon_1],[-10,15],'k')
grid on
% set(gca,'XLim', ax.XLim);
% set(gca,'YLim', ax.YLim);
hold off
xlabel('$\epsilon$','Interpreter','Latex');
ylabel('$V(\epsilon)$','Interpreter','Latex')
legend([h3,h2,h1],{'N = 10 ';['N = ',num2str(length(ind))];'N = \infty'},'Location','NorthWest');
axis([-20 1400 -10 30])
% print(gcf, '-depsc2', ['/Users/danielhurowitz/PROJ/NEG/Figs/gRing_ES_vs_cont.eps'])
%%
SpectDet=2*(cos(q) + 1/g * (q.^2+S^2/4)./(2*q) .*sin(q))-2;
figure;
hold on
plot(epsilon,SpectDet./((abs(epsilon_1-epsilon))/g))
plot([epsilon_1,epsilon_1],[-10,15],'k')
grid on