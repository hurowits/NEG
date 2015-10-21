%% Complexity saturation 1 (Fig 7a of NEG)
colors = 'rb'
figure;
axes('FontSize',24);
hold on;
grid on
h1=plot(linspace(0,10,100),squeeze(numComplex(:,1,:)),'r','LineWidth',2);
% h2=plot(linspace(0,5,100),squeeze(numComplex(:,2,:)),'b','LineWidth',2);
for iSigma = 1:2
    sigma = sigma_vec(iSigma)
    s_1_2=log(sinh(sigma/2)/(sigma/2))*2;
    s_1=log(sinh(sigma)/sigma);
    s_2=log(sinh(2*sigma)/(2*sigma))/2;
    
    a=exp((s-sigma)/2);
    b=exp((s+sigma)/2);
    max_x=fsolve(@(x)real(-log(-a+  x).* log(a./x) + log(b - x).* log(b./x)...
        - polylogN(2, 1 - a./x) + polylogN(2, 1 - b./x))/sigma-s/2,b/2);
    pcntComplex=log(max_x/a)/log(b/a);
    
    
    
    %     plot(s_vec(s_vec<s_1_2),zeros(size(s_vec(s_vec<s_1_2))),'k','LineWidth',3);
    %     plot(s_vec(s_vec>s_1_2),N*ones(size(s_vec(s_vec>s_1_2))),'k','LineWidth',3);
    %     plot(s_vec,N*pcntComplex*ones(size(s_vec)),'-.k','LineWidth',4)
    plot(s_vec,N*pcntComplex*ones(size(s_vec)),['k','--'],'LineWidth',4)
    %     plot([s_1_2 s_1_2],[0 N],'k','LineWidth',3);
    %     plot([s_1_2,s_1_2],[0, N],[colors(iSigma),'--'],'LineWidth',2);
    plot([s_1,s_1],[0, N],['k--'],'LineWidth',2);
    %     plot([s_1,s_1],[0, N],[colors(iSigma),'--'],'LineWidth',2);
    %     plot([sigma,sigma],[0, 1.2*N],[colors(iSigma),'--'],'LineWidth',2);
    
end
xlabel('s');
ylabel('# complex eigenvalues')

% print(gcf, '-depsc2', ['/Users/danielhurowitz/PROJ/NEG/Figs/numComplex_100.eps'])

%% Complexity saturation 2 (Fig 7b of NEG)
figure;
axes('FontSize',24);
hold on;
grid on
p1=plot(s_vec,squeeze(numComplex(1:30,1,:)),'r','LineWidth',2);
hold on;
% p2=plot(linspace(0,10,1e2),squeeze(numComplex(1:20,2,:)),'b','LineWidth',2);
xlabel('s');
ylabel('# complex eigenvalues')
legend([p1(1),p2(1)],{['\sigma=',num2str(sigma_vec(1))];['\sigma=',num2str(sigma_vec(2))]})
axis([0 5 0 100])
for iSigma = 1:2
    sigma = sigma_vec(iSigma)
    s_1_2=log(sinh(sigma/2)/(sigma/2))*2;
    s_1=log(sinh(sigma)/sigma);
    s_2=log(sinh(2*sigma)/(2*sigma))/2;
    
    a=exp((s-sigma)/2);
    b=exp((s+sigma)/2);
    max_x=fsolve(@(x)real(-log(-a+  x).* log(a./x) + log(b - x).* log(b./x)...
        - polylogN(2, 1 - a./x) + polylogN(2, 1 - b./x))/sigma-s/2,b*0.9);
    pcntComplex(iSigma)=log(max_x/a)/log(b/a);
    
    
    
    %     plot(s_vec(s_vec<s_1_2),zeros(size(s_vec(s_vec<s_1_2))),'k','LineWidth',3);
    %     plot(s_vec(s_vec>s_1_2),N*ones(size(s_vec(s_vec>s_1_2))),'k','LineWidth',3);
    %     plot(s_vec,N*pcntComplex*ones(size(s_vec)),'-.k','LineWidth',4)
    plot(s_vec,N*pcntComplex(iSigma)*ones(size(s_vec)),['k','--'],'LineWidth',4)
    %     plot([s_1_2 s_1_2],[0 N],'k','LineWidth',3);
    %     plot([s_1_2,s_1_2],[0, N],[colors(iSigma),'--'],'LineWidth',2);
    
    %     plot([s_1,s_1],[0, N],[colors(iSigma),'--'],'LineWidth',2);
    %     plot([sigma,sigma],[0, 1.2*N],[colors(iSigma),'--'],'LineWidth',2);
    
end
% print(gcf, '-depsc2', ['/Users/danielhurowitz/PROJ/NEG/Figs/numComplex_100_alpha.eps'])
fsolve(@(x)real(-log(-a+  x).* log(a./x) + log(b - x).* log(b./x)...
        - polylogN(2, 1 - a./x) + polylogN(2, 1 - b./x))/sigma-s/2+log(geomean(g)),b*0.9)
%% mu(s)
mu_vec=findMu(linspace(0,10,100),sigma);
figure;
axes('FontSize',24);
hold on;

plot(mu_vec,linspace(0,10,100),'LineWidth',4)
plot(mu_vec,sigma*ones(1,100),'--g','LineWidth',2);
plot(mu_vec,sigma^2/6*mu_vec,'--k','LineWidth',2)
xlabel('\mu');
ylabel('s_{\mu}');
legend(['s_{\mu}         ';...
    's_{\infty}      ';...
    's= \mu\sigma^2/6'],'Location','SouthEast')
axis([0 20 0 6])
% print(gcf, '-depsc2', ['/Users/danielhurowitz/PROJ/NEG/Figs/s_mu.eps'])
%% Generate gallery of plots of N(epsilon)
flagLegend=1;
for iS=[10,20,22,27,32,35,37,40,70,100]
    
    s=s_vec(iS);
    
    options = optimset('TolX',1e-11,'TolFun',1e-11);
    mu = fsolve(@(mu)s-1/mu*log(sinh(mu*sigma)/(mu*sigma)),10,options)
    
    epsilon_j=epsilon_j_mat(iS,:);
    x0=logspace(log10(min(abs(epsilon_j))),log10(max(epsilon_j)),1e4);
    z=(2304/sigma^4*x0);
    x=sort(epsilon_j);
    figure;
    axes('FontSize',24);
    hold on;
    grid on
    plot(sort((epsilon_j)),(1:N),'o');
    
    plot(x0,log(x0/exp((s-sigma)/2))*N/sigma,'--r','LineWidth',4)
    plot(x0,sqrt(x0)*N/pi,'--g','LineWidth',4);
    NumStates = N*sigma^2/(24*pi^2)*1./(besselj(mu,sqrt(z)).^2+bessely(mu,sqrt(z)).^2);
    plot(x0,NumStates,'-.k','LineWidth',4)
    plot([exp((s-sigma)/2) exp((s-sigma)/2)], [0.1 N],':k','LineWidth',4)
    xlabel('\epsilon');
    ylabel('N(\epsilon)');
    if(s<5)
        title(['s=',num2str(s),', \sigma=',num2str(sigma),', \mu=',num2str(mu)]);
    else
        title(['s=',num2str(s),', \sigma=',num2str(sigma),', \mu=\infty']);
    end
    if(flagLegend==1)
        legend(['numerics     ';...
            'log-box      ';...
            'free particle';...
            'J_{\mu}      '],'Location','NorthWest')
        flagLegend=0;
    end
    
    axtype(3);
    axis([x(1) 1e2 1 1e4])
    % print(gcf, '-depsc2', ['/Users/danielhurowitz/PROJ/NEG/Figs/NumStates_',num2str(iS),'.eps'])
    
end