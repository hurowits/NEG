%% Initialization

clear all
colors='rgbckymrgbckym'

options = optimset('TolX',1e-11,'TolFun',1e-11);
% N=500;
N_vec = 100:100:2000;[100,500,1000,1500,2000];
N_vec=100 ;
numRealizations=30;
N=N_vec(1);
randInd1 = randperm(N);
randInd2 = randperm(N);
% iRR=10;

sigma_vec=[3 5];
sigma_vec=2;
% sigma=0;
% alpha_vec =.7;0.5;[0.1, 0.5, 1, 2, Inf ]; [100,1.1,1,0.9]
alpha_vec=Inf;0.5;.39;.9;.5;.7;
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
            s_vec=[s_1_2*0.7,s_1_2,.9*s_1,s_2,2*sigma];
            %             s_vec=[linspace(0,1,10),linspace(1,5,5),20];;s_1*0.9
            %             s_vec=30;
            s_vec=linspace(0,10,200);
            %             s_vec=10
            % s_vec=0.9*s_1
            %             s_vec= 2*sigma;
            %             s_vec=s_1_2/2;[s_1_2/2,s_1,2*sigma];[0,s_1_2,s_1,s_2,2*sigma];[1,3,5];
            %             s_vec=2*sigma;
            %             s_vec = linspace(0,10,100);
            %             s_vec=2*sigma;
            %             s_vec=s_2;s_1_2/2
            %             s_vec=linspace(0,1e1,100);
            s_vec=5;
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
% a=exp((s-sigma)/2);
% b=exp((s+sigma)/2);
% max_x=fsolve(@(x)real(-log(-a+  x).* log(a./x) + log(b - x).* log(b./x)...
%     - polylogN(2, 1 - a./x) + polylogN(2, 1 - b./x))/sigma-s/2,b/2);
% pcntComplex=max_x/max(abs(epsilon_s));

% axis([s_vec(1) s_vec(end)  -max(real(gap)) -min(real(gap))/5])
% Inverse localization length
% nu=gap2(iS);
% iS=6;%real
% iS=21;%complex
% iS=27;%s_c
% iS=40;%s>s_c
% iS=50;

%%
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

%% Plot integrated density of states, (NEG, Fig. 2)
%     for iS=1:length(s_vec)
for iS=1:length(s_vec)
    s=s_vec(iS);
    x= epsilon_j_mat(iS,:);
    if (sigma==0)
        mu=0;
    else
        mu = fsolve(@(mu)(s)-1/mu*log(sinh(mu*sigma)/(mu*sigma)),1,options);
    end
    figure;
    set(gca,'FontSize',24)
    hold on;
    grid on
    plot(sort(epsilon_j_mat(iS,:)),(1:N)/N,'bo','LineWidth',4)
    
    %         if (s<sigma)
    if (mu==0)
        plot(x,1/N*log(x/x(1)),'--k','LineWidth',4)
    else
        plot(x,1/N*(x/x(1)).^mu,'--k','LineWidth',4)
    end
    %         else
    if (alpha>=1)
        plot(x,log(x/exp((s-sigma)/2))/sigma,'--r','LineWidth',4)
    else
        %             plot(x,log(x/exp((s-sigma)/2)/min(g))/((alpha+1)*sigma),'--r','LineWidth',2)
        %             plot(x,log(exp(1)*x/x(end)),'--r','LineWidth',4)
    end
    %         end
    plot(x,(x/x(end)).^mu_alpha,'--g','LineWidth',4)
    if (s==0)
        %             plot(x,1./log(x).^2,['r',':'])
        %           plot(x,1./log(x).^2,[colors(iC),':'])
    end
    iC = iC+1;
    xlabel('\epsilon')
    ylabel('N(\epsilon)')
    if (iS==1 && iAlpha==5)
        legend(['Numerics    ';...
            '\mu_s       ';...
            %                 'log-box     ';...
            '\mu_{\alpha}'],'Location','SouthEast');
    end
    axtype(3);
    %               title(['\alpha = ',num2str(alpha),', \mu_{\alpha} = ',num2str(mu_alpha),', \mu_s = ',num2str(mu)]);
    plot([exp((s-sigma)/2) exp((s-sigma)/2)], [1/N 1],':k','LineWidth',4)
    
    axis([x(1) x(end) 1/N 1])
    %      print(gcf, '-depsc2', ['/Users/danielhurowitz/PROJ/NEG/Figs/N_E_',num2str(iS),'_',num2str(iAlpha),'.eps'])
    %      print(gcf, '-depsc2', ['/Users/danielhurowitz/PROJ/NEG/Figs/N_E_0.eps'])
    %     end
    %% Plot spectral density
    %         iS=2
    %         s=s_vec(iS);
    %
    %         E_vec3=epsilon_j_mat(iS,:);
    %         figure;
    %
    %         axes('FontSize',24);
    %
    %         hold on
    %         grid on
    %         plot(0,0,'b','LineWidth',2)
    %         plot(0,0,'r','LineWidth',2)
    %         plot(0,0,'g','LineWidth',2)
    %         plot(0,0,'k','LineWidth',2)
    %         plot(sort(epsilon_j_mat(iS,:)),(1:N)/N,'.','LineWidth',2)
    %         mu = fsolve(@(mu)(s)-1/mu*log(sinh(mu*sigma)/(mu*sigma)),1,options);
    %
    %             plot(E_vec3,log(E_vec3/exp((s-sigma)/2))/sigma,':r','LineWidth',2)
    %             plot(E_vec3,log(E_vec3/E_vec3(1))/log(E_vec3(end)/E_vec3(1)),':r','LineWidth',2)
    %
    %
    % %         plot(E_vec3, (E_vec3/E_vec3(end)).^(mu),'r:','LineWidth',2)
    %         plot(E_vec3, 1/N*(E_vec3/E_vec3(1)).^(mu),'r:','LineWidth',2)
    %
    %         %         plot(E_vec3,log(E_vec3/exp((s-sigma)/2))/sigma,[colors(iAlpha),'.-'],'LineWidth',1)
    %
    % %         plot(E_vec3, 1/N*(E_vec3/E_vec3(1)).^(mu_alpha),'--g','LineWidth',2)
    %         plot(E_vec3, (E_vec3/E_vec3(N)).^(mu_alpha),'--g','LineWidth',2)
    %         z=(2304/sigma^4*E_vec3);
    % %         z=E_vec3;
    %         NumStates = sigma^2/(24*pi^2)*1./(besselj(mu,sqrt(z)).^2+bessely(mu,sqrt(z)).^2);
    %
    %         plot(z,NumStates,'--k','LineWidth',2)
    %         xlabel('\epsilon');
    %         ylabel('N(\epsilon)');
    %         legend(['numerics      ';...
    %             '\mu_s         ';...
    %             '{\mu_{\alpha}}';...
    %             'J_{\mu}       '],'Location','SouthEast')
    %         title(['\alpha=',num2str(alpha),', s=',num2str(s_vec(iS))]);
    %         axis([epsilon_j_mat(iS,1) epsilon_j_mat(iS,end) 1/N 1])
    % %                 legend(num2str(alpha_vec'))
    %         axtype(3)
end
%end %iAlpha
%% Plot electrostatic potential along real axis 
for iS=1
    
    s=s_vec(iS);
   
    [kappa,signKappa,expKappa,E_vec] = thoulessKappa(epsilon_j_mat(iS,1:end),N,g);
    kappa_p = kappa(signKappa>0);
    E_vec_p = E_vec(signKappa>0);
    [peaks_p,ind_p]=findpeaks(kappa_p);
    kappa_m = kappa(signKappa<0);
    E_vec_m = E_vec(signKappa<0);
    [peaks_m,ind_m]=findpeaks(kappa_m);
    %%
    figure;
    axes('FontSize',24);
    hold on;
    grid on
    
    %             plot(E_vec,kappa/N-sum(log(g))/N);
    plot(E_vec,signKappa.*(kappa/N))
    % plot(E_vec,(kappa/N),'r')
    %     plot([E_vec(signKappa>0)],[log(expKappa(signKappa>0))/N]-sum(log(g))/N,'b','LineWidth',2);
    %         h1=  plot(E_vec,kappa/N,':r','LineWidth',4);
    
    h1=  plot(E_vec(signKappa>0),kappa(signKappa>0)/N,':b','LineWidth',3);
    %     h1=  plot(E_vec(signKappa<0),kappa(signKappa<0)/N,':r','LineWidth',2);
    %         h1=  plot(E_vec(signKappa<0),kappa(signKappa<0)/N,'r','LineWidth',4);
    
    plot(E_vec_p(ind_p),(peaks_p)/N,'g','LineWidth',4)
    %     plot(E_vec_m(ind_m),-(peaks_m)/N,'--r','LineWidth',2)
    %             plot(E_vec,spline(E_vec(ind_p),peaks/N-sum(log(g))/N,E_vec),'--b','LineWidth',2)
    h2 = plot(E_vec,ones(1,length(E_vec))*s/2,'--r','LineWidth',4);
    %     plot(epsilon_j_mat(iS,:),ones(1,N)*s/2,'k*')
    h2=          plot(E_vec,ones(1,length(E_vec))*log(2*cosh(N*s/2)-2)/N,'--r','LineWidth',4);
    %     plot(epsilon_j_mat(iS,:),ones(1,N)*log(2*cosh(N*s/2)-2)/N,'k*')
    %           h3 =   plot(epsilon_j_mat(iS,:),zeros(1,N),'.k','LineWidth',.01);
    %         title(['s = ',num2str(s),', \mu = ',num2str(mu_vec(iS)), ', \sigma = ',num2str(sigma),', \alpha = ',num2str(alpha)]);
    %     plot(E_vec,(kappa/N),'r--')
    
    xlabel('\epsilon');
    ylabel('V(\epsilon)');
    
    legend([h1,h2],...
        {'V(\epsilon)';...
        'V(0)'},...
        'Location','SouthEast')
    hold off
end
%% Plot spectrum in complex plane
figure;
h=plot(-epsilon_s_mat(1,:),'.','MarkerSize',20);
% axis([-1e-2 200 -40 40])
set(gca,'FontSize',24);
% hold on; grid on;
xlabel('Re[\lambda]');ylabel('Im[\lambda]');

%% Scatter diagram of the spectrum
iS_max=100;
E=epsilon_s_mat(1:10:iS_max,:).';
C=round(repmat(s_vec(1:10:iS_max)',1,N)');
figure;scatter(-real(E(:)),imag(E(:)),30,(C(:)))
%%
for iS=1:length(s_vec)
    figure; plot(-epsilon_s_mat(iS,:),'.')
end
%% Analytical formula for V(epsilon) for Large S (weak disorder)
a=exp((s-sigma)/2);
b=exp((s+sigma)/2);
x=[linspace(eps,max(abs(epsilon_s_mat(iS,:))),5*1e2)];
kappa_large_S_1 =-log(abs(-a+x)).* log(a./x) + log(abs(b - x)).* log(b./x)...
    - polylogN(2, 1 - a./x) + polylogN(2, 1 - b./x);
%
% kappa_large_S_2 =-(pi^2/6) + log(b - x).*log(b./x) + log(x).* log(x./a) +...
%     polylogN(2, 1 - b./x) + polylogN(2, a./x);
%
% kappa_large_S_3 =log(b/a).*log(x) + polylogN(2, a./x) - polylogN(2, b./x);
%
% max_x=fsolve(@(x)real(-log(-a+  x).* log(a./x) + log(b - x).* log(b./x)...
%     - polylogN(2, 1 - a./x) + polylogN(2, 1 - b./x))/sigma-s/2,b/2)
% pcntComplex(iS)=max_x/b;
% %
%
plot(x,(kappa_large_S_1)/sigma,'--b','LineWidth',2)
% plot(x,real(kappa_large_S_2)/sigma,'--g','LineWidth',2)
% plot(x,real(kappa_large_S_3)/sigma,'--b','LineWidth',2)
%% Small S (strong disorder), "French" expressions
figure;
D0=sqrt(2*(cosh(s)-1)/s^2);
x0=logspace(-10,log10(max(E_vec)),1e5);
x=sqrt(384*D0.^3/sigma^4*x0);

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
% title(['s=',num2str(s),' \sigma=',num2str(sigma)]);
% title(['s=',num2str(s)]);

%
% figure(2);
% hold on
% %     plot(E_vec,real(log(expKappa))/N,colors(iSigma))
% plot([E_vec(signKappa>0)],[kappa(signKappa>0)/N],colors(iSigma),'LineWidth',2);
% %     plot(E_vec,ones(1,length(E_vec))*s/2+sum(log(g))/N,'--r','LineWidth',2);
%
% legend(num2str(sigma_vec'))
% end
% print(gcf, '-depsc2', ['/Users/danielhurowitz/PROJ/NEG/Figs/envelope_',num2str(iS),'.eps'])
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
%% CDF of Hermitian eigenvalues
iS=1;
% for iS=[1,10,20,22,27,32,35,37,40,70]
s=s_vec(iS);
Epsilon_n = 2*cosh(s/2)-2*cos((0:N-1)*2*pi/N);

options = optimset('TolX',1e-11,'TolFun',1e-11);
% mu = fsolve(@(mu)s-1/mu*log(sinh(mu*sigma)/(mu*sigma)),100,options)
mu=0
epsilon_j=epsilon_j_mat(iS,:);
x=linspace(min(epsilon_j),max(epsilon_j),10000);
x0=logspace(log10(min(abs(epsilon_j))),log10(max(epsilon_j)),1e4);

figure;
axes('FontSize',24);
hold on;
plot(sort((epsilon_j)),(1:N),'.-');
z=(2304/sigma^4*x0);
x=sort(epsilon_j);
plot(x0,log(x0/exp((s-sigma)/2))*N/sigma,'--r','LineWidth',2)
plot(x0,sqrt(x0)*N/pi,'--g','LineWidth',2);
plot(sort(Epsilon_n),1:N,'--g','LineWidth',2);

NumStates = N*sigma^2/(24*pi^2)*1./(besselj(mu,sqrt(z)).^2+bessely(mu,sqrt(z)).^2);
plot(x0,NumStates,'-.k','LineWidth',2)
%     plot(x0,x0.^mu,':c','LineWidth',2);
% plot(x0,N*(x0/x0(end)).^(1/(1+2*Delta)),':c','LineWidth',2);
% plot(x0,(x0/x0(1)).^(mu),':y','LineWidth',2);

%  plot(x,exp(-1./(sqrt(x)*sigma^2)),'r')
plot([exp((s-sigma)/2) exp((s-sigma)/2)], [0.1 N],':k','LineWidth',2)
xlabel('\epsilon');
ylabel('N(\epsilon)');
title(['s=',num2str(s),', \sigma=',num2str(sigma),' \mu=',num2str(mu)]);
legend(['numerics     ';...
    'log-box      ';...
    'free particle';...
    'J_{\mu}      ']);
axtype(3);
%     axis([x(1) x(end) 1 1e4])
% end
% print(gcf, '-depsc2', '/Users/danielhurowitz/PROJ/NEG/Figs/CDF_eps_weak.eps')
% print(gcf, '-depsc2', '/Users/danielhurowitz/PROJ/NEG/Figs/CDF_log_eps_weak.eps')
%%
%% CDF of Hermitian eigenvalues with off diagonal disorder
iS=2;
for iS=1:length(s_vec)
    s=s_vec(iS);
    Epsilon_n = 2*cosh(s/2)-2*cos((0:N-1)*2*pi/N);
    
    options = optimset('TolX',1e-11,'TolFun',1e-11);
    mu = fsolve(@(mu)s-1/mu*log(sinh(mu*sigma)/(mu*sigma)),100,options)
    %     mu=0
    epsilon_j=epsilon_j_mat(iS,:);
    x=linspace(min(epsilon_j),max(epsilon_j),10000);
    x0=logspace(log10(min(abs(epsilon_j))),log10(max(epsilon_j)),1e4);
    
    figure;
    axes('FontSize',24);
    hold on;
    plot(sort((epsilon_j)),(1:N),'.-');
    z=(2304/sigma^4*x0);
    x=sort(epsilon_j);
    %       plot(x0,x0.^mu,':c','LineWidth',2);
    plot(x0,N*(x0/x0(end)).^(alpha/(alpha+1)),':g','LineWidth',2);
    plot(x0,4*(x0/x(4)).^(mu),':r','LineWidth',2);
    plot(sort(Epsilon_n),1:N)
    %  plot(x,exp(-1./(sqrt(x)*sigma^2)),'r')
    %     plot([exp((s-sigma)/2) exp((s-sigma)/2)], [0.1 N],':k','LineWidth',2)
    xlabel('\epsilon');
    ylabel('N(\epsilon)');
    title(['s=',num2str(s),', \sigma=',num2str(sigma),' \mu=',num2str(mu)]);
    legend(['numerics     ';...
        '\alpha       ';...
        '\mu          '])
    axtype(3);
end
%     axis([x(1) x(end) 1 1e4])
% print(gcf, '-depsc2', '/Users/danielhurowitz/PROJ/NEG/Figs/CDF_eps_weak.eps')
% print(gcf, '-depsc2', '/Users/danielhurowitz/PROJ/NEG/Figs/CDF_log_eps_weak.eps')