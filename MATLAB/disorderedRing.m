%% Initialization

clear all
colors='rgbckymrgbckym'

options = optimset('TolX',1e-11,'TolFun',1e-11);
% N=500;
N_vec = 100:100:2000;[100,500,1000,1500,2000];
N_vec=100;
numRealizations=1;
N=N_vec(1);
randInd1 = randperm(N);
randInd2 = randperm(N);
% iRR=10;

% sigma_vec=[3 5];
sigma_vec=0;
% sigma=0;
% alpha_vec =.7;0.5;[0.1, 0.5, 1, 2, Inf ]; [100,1.1,1,0.9]
alpha_vec=1.5;0.3;0.5;0.5;.39;.9;.5;.7;
% alpha_vec = 0.5;
% alpha_vec=2
% alpha_vec=0.4
% alpha_vec=0.8
iC=1;
%%
% for iAlpha=1:length(alpha_vec)
iAlpha=1;
for iSigma=1:length(sigma_vec)
    for iR = 1:numRealizations
        randInd1 = randperm(N);
        randInd2 = randperm(N);
        iR
        for iN=1:length(N_vec)
            N=N_vec(iN);
            display(['N = ',num2str(N)])
            x0 = linspace(-1,1,N);
            % x0 = 2*linspace(0,1,N) - 1;
            %         B0 = linspace(-1,0,N);
            %         randInd1 = randperm(N);
            %         randInd2 = randperm(N);
            
            x0 = x0(randInd1) - mean(x0);
            x_n = (1:N)/N;
            %             x_n = rand(1,N)
            x_n = x_n(randInd2);
            
            sigma=sigma_vec(iSigma);
            
            alpha=alpha_vec(iAlpha);
            g=(x_n./N).^(1/alpha);
            
            %             g(1)=1e-2;0.9;
            % %             g(5)=1e-2
            %             g(6)=1e-2;
            %             g(12)=1e-2;
            %                         g=ones(1,N);
            if(alpha>1)
                mu_alpha=1/2;
            else
                mu_alpha=alpha/(alpha+1);
            end
            s_1_2=log(sinh(sigma/2)/(sigma/2))*2;
            s_1=log(sinh(sigma)/sigma);
            s_2=log(sinh(2*sigma)/(2*sigma))/2;
            
            %             hwb=waitbar(0,'Please wait');
            %             mu_vec=[0.3,0.5,0.7,1,1.5,2,4]
            %             mu_vec=200;
            %             s_vec=1./mu_vec.*log(sinh(sigma*mu_vec)./sigma./mu_vec);
            s_vec=[s_1_2*0.7,.9*s_1,2*sigma];
%             s_vec=s_1_2*0.5
            %             s_vec=[linspace(0,1,10),linspace(1,5,5),20];;s_1*0.9
%             %             s_vec=30;
%             s_vec=linspace(0,10,200);
            %             s_vec=10
%             s_vec=0.8*s_1
            %             s_vec= 2*sigma;
            %             s_vec=s_1_2/2;[s_1_2/2,s_1,2*sigma];[0,s_1_2,s_1,s_2,2*sigma];[1,3,5];
            %             s_vec=2*sigma;
            %             s_vec = linspace(0,10,100);
            %             s_vec=2*sigma;
            %             s_vec=s_2;s_1_2/2
            %             s_vec=linspace(0,1e1,100);
            s_vec=linspace(0,10/N,1e2);
            for iS=1:length(s_vec)
                %                 mu=mu_vec(iS);
                
                
                s = s_vec(iS);
                %     s=4;
                x=x0*sigma+s;
                
                
                [epsilon_j,gamma] = SolveNHR2(x,g); %find eigenvalues without off diagonal asymmetry
                epsilon_j_mat(iS,:)=epsilon_j;
                gamma_mat(iS,:)=gamma;
                %                 [epsilon_s,Drift(iS),Diffusion(iS)] = SolveNHR(x,g);  %find true eigenvalues
                [epsilon_s] = SolveNHR(x,g);  %find true eigenvalues
                epsilon_s_mat(iS,:)=epsilon_s;
                numComplex(iR,iSigma,iS)=sum(imag(epsilon_s)~=0); %how many complex eigenvalues?
                
                gap(iS)=epsilon_s(N-1);
                
            end
            
            %             close(hwb);
            
        end
    end
end

%% Plot integrated density of states for sigma=0

figure;
axes('FontSize',24)
hold on

h1=plot(sort(epsilon_j_mat),(1:N)/N,'*');
h2=plot(sort(epsilon_j_mat(1,:)),(epsilon_j_mat(1,:)/epsilon_j_mat(1,end)).^mu_alpha,'--b','LineWidth',2)
% h3=plot(sort(epsilon_j_mat(end,:)),(epsilon_j_mat(end,:)/epsilon_j_mat(end,end)).^0.5,'--r','LineWidth',2)
h3=plot(sort(epsilon_j_mat(end-1,:)),1/N*(epsilon_j_mat(end-1,:)/epsilon_j_mat(end-1,1)).^alpha,'--r','LineWidth',2)
h3=plot(sort(epsilon_j_mat(2,:)),1/N*(epsilon_j_mat(2,:)/epsilon_j_mat(2,1)).^alpha,'--r','LineWidth',2)
axtype(3)
axis([1e-15 1e-5 1e-4 1])
legend([h2,h3],...
    {['\mu_{\alpha}=',num2str(mu_alpha)];...
    '\mu=\alpha'},...
    'Location','SouthEast')

%     for iS=1:length(s_vec)
%% Plot N(epsilon), Fig. 2 of NEG
for iS=1:length(s_vec)
    s=s_vec(iS);
    x= epsilon_j_mat(iS,:);
    if (sigma==0)
        mu=0;
    else
        mu = fsolve(@(mu)(s)-1/mu*log(sinh(mu*sigma)/(mu*sigma)),1,options);
    end
    
    z = 144*16*x/sigma^4;
    NumStates = N*sigma^2/(24*pi^2)./(besselj(mu,sqrt(z)).^2+bessely(mu,sqrt(z)).^2);
    
    figure;
    set(gca,'FontSize',24)
    hold on;
    grid on
    h1=plot(sort(epsilon_j_mat(iS,:)),(1:N)/N,'bo','LineWidth',4)
        h5=plot(x(1:end),NumStates(1:end)./(NumStates(1))*1/N,'-r','LineWidth',2)

    %         if (s<sigma)
    if (mu==0)
       h2= plot(x,1/N*log(x/x(1)),'--k','LineWidth',4)
    else
        h2=plot(x(1:end),1/N*(x(1:end)/x(1)).^mu,'--k','LineWidth',4)
    end
    %         else
    if (alpha>=1)
        h4=plot(x,log(x/exp((s-sigma)/2))/sigma,'--.r','LineWidth',4)
    else
        %             plot(x,log(x/exp((s-sigma)/2)/min(g))/((alpha+1)*sigma),'--r','LineWidth',2)
        %             plot(x,log(exp(1)*x/x(end)),'--r','LineWidth',4)
    end
    %         end
   h3= plot(x,(x/x(end)).^mu_alpha,'--g','LineWidth',4)
    if (s==0)
        %             plot(x,1./log(x).^2,['r',':'])
        %           plot(x,1./log(x).^2,[colors(iC),':'])
    end
    plot([exp((s-sigma)/2) exp((s-sigma)/2)], [1/N 1],':k','LineWidth',4)
    %
    iC = iC+1;
    xlabel('$\epsilon$','Interpreter','Latex')
    ylabel('$N(\epsilon)$','Interpreter','Latex')
    if (iS==1 && iAlpha==1)
%         legend_h=legend([h1,h2,h3,h5,h4],...
%         {'Numerics';...
%          '$\epsilon^{\mu_s}$';...
%          '$\sqrt{\epsilon}$';...
%          '$J_{\mu_s}(\epsilon)$';...
%          'log-box'},'Location','NorthWest');
%      set(legend,'Interpreter','Latex')

legend_h=legend([h1,h2,h3,h5],...
        {'Numerics';...
         '$\epsilon^{\mu_s}$';...
         '$\epsilon^{\mu_{\alpha}}$';...
         '$J_{\mu_s}(\epsilon)$'},'Location','SouthEast');
     set(legend,'Interpreter','Latex')
%      set(legend_h,'Color','none')
    end
    axtype(3);
    title(['$\mu_s$ = ',num2str(mu)],'Interpreter','Latex');
    axis([x(1) x(end) 1/N 1])
%          print(gcf, '-depsc2', ['/Users/danielhurowitz/PROJ/NEG/Figs/N_E_',num2str(iS),'_',num2str(iAlpha),'.eps'])
% %          print(gcf, '-depsc2', ['/Users/danielhurowitz/PROJ/NEG/Figs/N_E_',num2str(iS),'_',num2str(iAlpha),'_french.eps'])
% print(gcf, '-depsc2', ['/Users/danielhurowitz/PROJ/NEG/Figs/N_E_',num2str(iS),'_french.eps'])
    %      print(gcf, '-depsc2', ['/Users/danielhurowitz/PROJ/NEG/Figs/N_E_0.eps'])
    %     end
    
end
% save('N_E_french.mat');
%end %iAlpha
%% Plot electrostatic potential along real axis
for iS=7:8%1:length(s_vec)
    
    s=s_vec(iS);
%     E_vec=linspace(epsilon_j_mat(iS,1),epsilon_j_mat(iS,end)*1.1,1e5);
E_vec=linspace(0,epsilon_j_mat(iS,10)*1.1,1e5);
    [kappa,signKappa,expKappa,E_vec] = thoulessKappa(epsilon_j_mat(iS,1:end),N,g,E_vec);
    kappa_p = kappa(signKappa>0);
    E_vec_p = E_vec(signKappa>0);
    [peaks_p,ind_p]=findpeaks(kappa_p);
    kappa_m = kappa(signKappa<0);
    E_vec_m = E_vec(signKappa<0);
    [peaks_m,ind_m]=findpeaks(kappa_m);
    %%
    a=1 %stretching factor
    figure;
    axes('FontSize',24);
    hold on;
    grid on
    
    %             plot(E_vec,kappa/N-sum(log(g))/N);
    plot(E_vec.^a,(kappa/N))
    % plot(E_vec,(kappa/N),'r')
    %     plot([E_vec(signKappa>0)],[log(expKappa(signKappa>0))/N]-sum(log(g))/N,'b','LineWidth',2);
    %         h1=  plot(E_vec,kappa/N,':r','LineWidth',4);
    
    h1=  plot(E_vec(signKappa>0).^a,kappa(signKappa>0)/N,'.b','LineWidth',3);
    %     h1=  plot(E_vec(signKappa<0),kappa(signKappa<0)/N,':r','LineWidth',2);
    %         h1=  plot(E_vec(signKappa<0),kappa(signKappa<0)/N,'r','LineWidth',4);
    
%     h5=plot(E_vec_p(ind_p).^a,(peaks_p)/N,'g','LineWidth',4)
%     set(h5,'Color',[0 1 0])
%     h5=plot(E_vec.^a,kappa/N,'g','LineWidth',4)
    %     plot(E_vec_m(ind_m),-(peaks_m)/N,'--r','LineWidth',2)
    %             plot(E_vec,spline(E_vec(ind_p),peaks/N-sum(log(g))/N,E_vec),'--b','LineWidth',2)
%     h2 = plot(E_vec.^a,ones(1,length(E_vec))*s/2,'--r','LineWidth',4);
    %     plot(epsilon_j_mat(iS,:),ones(1,N)*s/2,'k*')
    h2=          plot(E_vec.^a,ones(1,length(E_vec))*log(2*cosh(N*s/2)-2)/N,'--r','LineWidth',4);
    %     plot(epsilon_j_mat(iS,:),ones(1,N)*log(2*cosh(N*s/2)-2)/N,'k*')
    %           h3 =   plot(epsilon_j_mat(iS,:),zeros(1,N),'.k','LineWidth',.01);
    %         title(['s = ',num2str(s),', \mu = ',num2str(mu_vec(iS)), ', \sigma = ',num2str(sigma),', \alpha = ',num2str(alpha)]);
    %     plot(E_vec,(kappa/N),'r--')
    set(h1,'Color',[0 0 1])
set(h2,'Color',[1 0 0])
% set(h5,'Color',[0 .7 0])
    xlabel('$\epsilon$','Interpreter','Latex');
    ylabel('$V(\epsilon)$','Interpreter','Latex');
    
%     legend([h5,h2],...
%         {'$V(\epsilon)$';...
%         '$V(0)$'},...
%         'Location','SouthEast','Interpreter','Latex')
    hold off
end
%% Plot spectrum in complex plane
figure;
h=plot(-epsilon_s_mat(2,:),'.','MarkerSize',20);
% axis([-1e-2 200 -40 40])
set(gca,'FontSize',24);
% hold on; grid on;
xlabel('Re[\lambda]');ylabel('Im[\lambda]');
axis([0 6 -0.12 0.12])
%% Scatter diagram of the spectrum. Color corresponds to s
iS_max=100;
E=epsilon_s_mat(1:10:iS_max,:).';
C=round(repmat(s_vec(1:10:iS_max)',1,N)');
figure;scatter(-real(E(:)),imag(E(:)),30,(C(:)))

%% Analytical formula for V(epsilon) for Large S (weak disorder)
a=exp((s-sigma)/2);
b=exp((s+sigma)/2);
x=[linspace(eps,max(abs(epsilon_s_mat(iS,:))),5*1e2)];
kappa_large_S_1 =-log(abs(-a+x)).* log(a./x) + log(abs(b - x)).* log(b./x)...
    - polylogN(2, 1 - a./x) + polylogN(2, 1 - b./x);
plot(x,(kappa_large_S_1)/sigma,'--b','LineWidth',2)
% plot(x,real(kappa_large_S_2)/sigma,'--g','LineWidth',2)
% plot(x,real(kappa_large_S_3)/sigma,'--b','LineWidth',2)
%% Small S (strong disorder), "French" expressions
figure;
D0=1;sqrt(2*(cosh(s)-1)/s^2);
x0=logspace(-10,log10(max(epsilon_j_mat(iS,:))),1e5);
x=sqrt(144*16*D0.^3/sigma^4*x0);

options = optimset('TolX',1e-11,'TolFun',1e-11);
mu = fsolve(@(mu)(s)-1/mu*log(sinh(mu*sigma)/(mu*sigma)),1,options);
% mu=2

kappa_small_S =  -0.5*x.* (besselj(mu,x).*(besselj(mu-1,x)-besselj(mu+1,x)) + ...
    bessely(mu,x).*(bessely(mu-1,x)-bessely(mu+1,x)))./...
    (besselj(mu,x).^2+bessely(mu,x).^2);


% plot(x0,kappa_small_S*sigma^2/(48*D0^2)*s,'--y','LineWidth',2)
% plot(x0,kappa_small_S/mu * s/2,'--y','LineWidth',2)
plot(x0,kappa_small_S/(2*mu)*s,'--c','LineWidth',2)

% a=get(gca);
% plot(E_vec,sigma^2/24./(4-(E_vec+2*sinh(sigma/2)/sigma).^2),'--r','LineWidth',2)

% plot(E_vec,4*sigma./(E_vec),'k')
% axis([a.XLim a.YLim])

xlabel('\lambda');
ylabel('\kappa');
%% Histogram of Hermitian eigenvalues
figure;
axes('FontSize',24);
hold on;
[h,x_h]=hist(((epsilon_j)),50);
bar(x_h,h/N);
bar(x_h(1),h(1)/N,'r')
dx=x_h(2)-x_h(1);

[h2,x_h2]=hist(epsilon_j,100);

x=linspace(x_h(1),x_h(end),1000);
plot(x,dx./(sigma*x),'r','LineWidth',2);
% plot(x,dx*exp(-pi./(sqrt(x)*sigma)).*pi/sigma.* x.^(-3/2)/2,'r','LineWidth',2);

xlabel('\epsilon_j');
ylabel('p(\epsilon_j)');
title(['s=',num2str(s),', \sigma=',num2str(sigma)]);
% print(gcf, '-depsc2', '/Users/danielhurowitz/PROJ/NEG/Figs/P_E_small_S.eps')
% print(gcf, '-depsc2', '/Users/danielhurowitz/PROJ/NEG/Figs/P_eps_strong.eps')
