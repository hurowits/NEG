%%
N=50
k=2*pi/N * (0:N-1);
R_vec = linspace(-1,5,1e4);
s=0;
for iR = 1:length(R_vec)
    R=R_vec(iR);
    V(iR) =(N/2/pi) * sum(log((2*cosh(s/2)-2*cos(k)-R)));
    V2(iR) = prod(2*cosh(s/2)-2*cos(k)-R);
end
R_vec2=linspace(-1,0,1e3);

E=linspace(2*cosh(s/2)-2,2*cosh(s/2)+2,100);
dE=E(2)-E(1)
rho = N/4/pi ./sqrt(1-(cosh(s/2)-E/2).^2);
% plot(E,rho*dE,E,N/4/pi* dE*(1./sqrt(E)))


figure;
axes('FontSize',24);
hold on
plot(E,dE*rho,'-k','LineWidth',4)
% plot(R_vec,real(V)/N,R_vec2,log(2)/pi + N/2/pi*sqrt(abs(R_vec2)),R_vec,ones(1,1e3)*log(2)/pi,':','LineWidth',4);
plot(R_vec,real(V2),'LineWidth',4);
xlabel('R');
% ylabel('V(R)/N');
axis([-.2 4.2 -5 5]);
%legend({'\rho (R)';'V(R)';' |R|^{1/2}'})
% legend({'\rho (R)';'V(R)'})
% axis([-.2 4.2 -2 2]);
% print(gcf, '-depsc2', ['/Users/danielhurowitz/PROJ/NEG/Figs/V_R_stepbystep_continuum.eps'])


%%
a=1;b=3;
r=linspace(-1,5,1e4);
figure;
axes('FontSize',24);
hold on
plot(linspace(a,b,10),ones(1,10)/N,'k','LineWidth',4)
plot(r,(r-a).*log(abs(r-a))-(r-b).*log(abs(r-b))+a-b,'LineWidth',4)
% legend({'\rho (R)';'V(R)'})
xlabel('R')
axis([-.2 4.2 -2.1 2]);
% print(gcf, '-depsc2', ['/Users/danielhurowitz/PROJ/NEG/Figs/V_R_stepbystep_uniform.eps'])

%%
r=linspace(-1,5,1e4);
figure;
axes('FontSize',24);
hold on
plot(2,1/N,'*k','LineWidth',4)
plot(r,log(abs(r-2)),'LineWidth',4)
legend({'\rho (R)';'V(R)'})
axis([-.2 4.2 -2.1 2]);xlabel('R')

% print(gcf, '-depsc2', ['/Users/danielhurowitz/PROJ/NEG/Figs/V_R_stepbystep_point.eps'])

xlabel('R')