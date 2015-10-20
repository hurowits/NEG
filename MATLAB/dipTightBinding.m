clear all
flag=1;
N=100;
n=0:N-1;



sigma=0;
x0=zeros(1,N);
% x0(N/2)=-sigma;
% x0(N/2+1)=sigma;

% x0(2)=-2*sigma
% x0(3)=2*sigma

% x0(20)=-sigma
% x0(21)=sigma
% x0(59)=sigma;
% x0(60)=-sigma;
g=ones(1,N);
% g(1)=1e-4;1e-1;0.5/N;1;N/100;0.00005;

    randInd2 = randperm(N);
alpha=Inf;
x_n = (1:N)/N;
    x_n = x_n(randInd2);
    g=(x_n./N).^(1/alpha);
% g(1)=1
g(1)=1;3*1e-4;
% g(10)=1e-9;
% g(15)=1e-9
% g(25)=1e-9
% g(11)=3*1e-5
% g(16)=3*1e-6
% g(15)=1e-3% g(70)=1e-4;
% g(90)=1e-8;
% x0=ones(1,N);
s_vec = linspace(0,1,100);
% s_vec=.01;
% s_vec=0.06;
% s_vec=[0.1,0.2]
% s_vec=0.35
for iS = 1:length(s_vec)
    
    s = s_vec(iS);
    x = x0 + s;
    %     g=exp(B0);
    %     g=ones(1,N);
    
    [epsilon_j,gamma,v] = SolveNHR2(x,g); %find eigenvalues without off diagonal asymmetry
    eigenStates(iS,:,:)=v;
    epsilon_j_mat(iS,:)=epsilon_j;
    gamma_mat(iS,:) = gamma;
    epsilon_s = SolveNHR(x,g);  %find true eigenvalues
    epsilon_s_mat(iS,:)=epsilon_s;
    
    %     lambda_s(iS,:) = SolveNHR(x,g);
    %     lambda(iS,:) =2*( cosh(s/2)*(cos(2*pi*n/N)-1)+1i*sinh(s/2)*sin(2*pi*n/N));
    S(iS)=sum(x)/N;
%     figure(1);
%     imagesc(epsilon_j,1:N,abs(v).^2,[0 0.04])
%     figure(1);
%     stem(1:N,abs(v(:,1)))
    %     E_vec(iS,:) = [linspace(0,0,100),linspace(0,max(abs(epsilon_s))*1.2,1e4)];
    %     %         E_vec = [linspace(-1,0,1e5),linspace(1e-5,1,1e5)];
    %     %             E_vec(iS,:) = logspace(-5,log(max(abs(epsilon_s))*1.2),1e5);
    % %
    %     for iE=1:length(E_vec)
    %         E=E_vec(iS,iE);
    %         kappa(iS,iE) = sum(log(abs(E-epsilon_j))-sum(log(g)));
    %         expKappa(iS,iE)=prod((E-epsilon_j));
    %     end
    
    %     figure(2);
    % %
    %     plot(epsilon_j,zeros(size(epsilon_j)),'*')
    %     axis([0 20 -1 1])
    %     title(['s=',num2str(s)])
    
    
end
%%
for iS=10
   
 s=s_vec(iS);
    %         % nu=imag(gap(iS));
    %     E_vec = [linspace(0,max(abs(epsilon_s_mat(iS,:))),1e4)]';
    %         E_vec = [linspace(0,epsilon_j_mat(iS,4),1e4)]';
    %     g=ones(1,N);
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
%     plot(E_vec,signKappa.*(kappa/N- sum(log(g))/N))
% plot(E_vec,(kappa/N),'r')
    %     plot([E_vec(signKappa>0)],[log(expKappa(signKappa>0))/N]-sum(log(g))/N,'b','LineWidth',2);
        h1=  plot(E_vec,kappa/N,':r','LineWidth',4);
    
    h1=  plot(E_vec(signKappa>0),kappa(signKappa>0)/N,'.b','LineWidth',4);
%     h1=  plot(E_vec(signKappa<0),kappa(signKappa<0)/N,':r','LineWidth',2);
%         h1=  plot(E_vec(signKappa<0),kappa(signKappa<0)/N,'r','LineWidth',4);
    
    plot(E_vec_p(ind_p),(peaks_p)/N,'--g','LineWidth',2)
%     plot(E_vec_m(ind_m),-(peaks_m)/N,'--r','LineWidth',2)
    %             plot(E_vec,spline(E_vec(ind_p),peaks/N-sum(log(g))/N,E_vec),'--b','LineWidth',2)
%     h2 = plot(E_vec,ones(1,length(E_vec))*s/2,'--r','LineWidth',4);
    h2=          plot(E_vec,ones(1,length(E_vec))*log(2*cosh(N*s/2)-2)/N,'--r','LineWidth',4);
    plot(epsilon_j_mat(iS,:),ones(1,N)*log(2*cosh(N*s/2)-2)/N,'k*')
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
%%
if (flag==1)
iS=1
s=s_vec(iS)

[kappa2,signKappa,expKappa,E2]=thoulessKappa(epsilon_j_mat(iS,1:end),N,g);
%%
figure;
axes('FontSize',24);
% set(h,'NextPlot','add')
hold on;
plot(E2,kappa2/N,'LineWidth',2);%+log(abs(E'-epsilon_j(1)))/N,'b')
E=linspace(epsilon_j(2),epsilon_j(end-1),100);
dE=E(2)-E(1)
rho = N/4/pi ./sqrt(1-(cosh(s/2)-E/2).^2);
plot(E,rho*dE,E,N/4/pi* dE*(1./sqrt(E)))
% 
% [kappa,signKappa,expKappa,E]=thoulessKappa(epsilon_j_mat(iS,2:end),100,g);
% plot(E,kappa/N)
% plot(E,log(abs(E'-epsilon_j(1)))/N,'c')
% plot(E,ones(1,length(E))*log(2*cosh(s*N/2)-2)/N,'r')
% plot(E,(sqrt(abs(E-epsilon_j(2)))),'--')
% plot(E,log(4*abs(E-epsilon_j(1))/min(g))/N+(sqrt(abs(E-epsilon_j(2)))),'--')
xlabel('\epsilon')
ylabel('V(\epsilon)')
% axis([min(E2) max(E2), -0.3 0.3]) 
end
% [kappa3,signKappa,E]=thoulessKappa(epsilon_j_mat(iS,[1,3:end]),3,g);
% plot(E,kappa3/N,'k')
% figure;plot(E,sqrt(abs(E-1e-5)),E,.3*sqrt(abs(E-1e-5)))
% figure;plot(E,sqrt(abs(E-epsilon_j(2)))-(kappa2'/N))
%%
[s,n]=meshgrid(s_vec,0:N-1);
lambda_cleanRing = 2*cosh(s/2).*(cos(n*2*pi/N)-1)+...
    2*1i*sinh(s/2).*sin(n*2*pi/N);
figure;
axes('FontSize',24);

hold on
plot(0,0,'b',0,0,'r',0,0,'k')
plot(s_vec,epsilon_j_mat,'-r','LineWidth',2);
% plot(s_vec,abs(real(epsilon_s_mat)),'b','LineWidth',2)
plot(s_vec,gamma_mat,'--b','LineWidth',2);

if (sigma~=0)
    
    
    plot(s_vec, exp(-sigma/2)*2*cosh(s_vec/2),'-c','LineWidth',2);
    plot(s_vec, exp(-s_vec/2)+exp(s_vec/2+sigma/2),'--c','LineWidth',2);
    plot(s_vec, exp(s_vec/2+sigma/2)+exp(-s_vec/2),'-c','LineWidth',2);
    plot(s_vec,2*cosh(s_vec/2)+2,'--m','LineWidth',2);
    plot(s_vec,2*cosh(s_vec/2)-2,'--m','LineWidth',2);
    % plot(s_vec,(2*cosh(s_vec/2)-2)exp(-sigma/2),'--m','LineWidth',2);
    
else
    plot(s_vec,min(g)*exp(s_vec/2)+exp(-s_vec/2),'--c','LineWidth',2);
    plot(s_vec,min(g)*exp(-s_vec/2)+exp(s_vec/2),'--c','LineWidth',2);
%     plot(s_vec,0.5*(2+min(g)-sqrt(4+min(g)^2))*ones(1,length(s_vec))/N/4);
%     plot(s_vec,min(g)*(2*cosh(s_vec/2)-2)*N)
    % plot(s_vec,2*cosh(s_vec/2)+2,'--m','LineWidth',2);
%     plot(s_vec,2*cosh(s_vec/2)-2,'--m','LineWidth',2);
    % plot(s_vec,N*min(g)*2*(cosh(s_vec/2)-1),'--m','LineWidth',2);
end
% plot(s_vec,-real(epsilon_s_mat),'--b','LineWidth',2);

% plot(s_vec,(s_vec.^2/4)./cosh(s_vec/2)/2/N,'g')
% plot(s_vec,(s_vec.^2/4+pi^2/N^2),'--g','LineWidth',2);%./cosh(s_vec/2),'--g')
% plot(s_vec,(s_vec.^2/4)*min(g)*N^2,'--c')
% plot(s_vec,min(g)*sinh(s_vec/2).^2./cosh(s_vec/2)*N,'y');
% plot(s_vec,N*(s_vec.^2/4).*min(g)/N,'--y','LineWidth',2)
% plot(s_vec,s_vec.^2/4/(sigma/2) ,'--b','LineWidth',2)
% plot(s_vec,-lambda_cleanRing','--g')
% plot(s_vec,2*cosh(s_vec/2),'--c')
% plot(s_vec,2*cosh(s_vec/2)+2,'--c','LineWidth',2);
% plot(s_vec,2*cosh(s_vec/2)-2,'--c','LineWidth',2);
% plot(s_vec,exp(s_vec/2)+exp(-(s_vec+sigma)/2),'--c','LineWidth',2);
% plot(s_vec,exp(-s_vec/2)+exp((s_vec-sigma)/2),'--c','LineWidth',2);
% plot(s_vec,2*cosh(s_vec/2)*exp((sigma)/2),'--m','LineWidth',2);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% plot(s_vec,exp(s_vec/2+sigma/2),'k')
% plot(s_vec,exp((s_vec-sigma)/2)/N^2,'k')
% plot(s_vec,N*s_vec.^2/4.*min(g)/(1-min(g)))
% plot(s_vec,N*(2*cosh(s_vec/2)-2).*min(g)/(1-min(g)));
plot(s_vec,N*s_vec.^2/4.*min(g)/(1-min(g)))
plot(s_vec,s_vec.*min(g)/(1-min(g)),'--g')
plot(s_vec,s_vec.^2/4+pi^2/N^2)
% plot([sigma/2 sigma/2],[ 0 15],'--k')
s_c = 2*acosh(1/(1-exp(-sigma/2)));
% s_c=2*N^(1/(M-1))
s_c=exp(-sigma/4)*2*sqrt(2)
% plot([s_c s_c],[ 0.001 15],'--k')
legend(['\gamma_n  ';...
    '\epsilon_n'])
   
    axtype(3)
%%
iS=1;
s=s_vec(iS)
D=cosh(s/2)
[kappa,signKappa,E]=thoulessKappa(epsilon_j_mat(iS,:),3,g);
E_max = epsilon_j_mat(iS,end);
figure;
hold on
% plot(E,(D*4*pi^2*sqrt(-(E-epsilon_j_mat(iS,2)))*pi-2),'b');

plot(E,(2*(D*4*pi^2)*sqrt(-(E-epsilon_j_mat(iS,2))/E_max)*pi-2)+log(abs(E-epsilon_j_mat(iS,1))/min(g)),'b');
plot(E,kappa,'g');
plot(E,log((2*cosh(N*s_vec(iS)/2)-2))*ones(1,length(E)),'--r')


iS=20;
s=s_vec(iS)
D=cosh(s/2)

[kappa,signKappa,E]=thoulessKappa(epsilon_j_mat(iS,:),10,g);
E_max = epsilon_j_mat(iS,end);
figure;
hold on
% plot(E,(D*4*pi^2*sqrt(-(E-epsilon_j_mat(iS,2)))*pi-2),'b');

plot(E,((D*4*pi^2)*sqrt(-(E-epsilon_j_mat(iS,2))/E_max)*pi-2)+log(abs(E-epsilon_j_mat(iS,1))/min(g)),'b');
plot(E,kappa,'g');
plot(E,log((2*cosh(N*s_vec(iS)/2)-2))*ones(1,length(E)),'--r')
%%
figure;
axes('FontSize',24);
hold on
for iN=1:N
    
    I1=imag(epsilon_s_mat(:,iN))==0;
    I2=imag(epsilon_s_mat(:,iN))~=0;
    hold on
    plot(S(I1),(-real(epsilon_s_mat(I1,iN))),'k-','LineWidth',2');
    plot(S(I2),(-real(epsilon_s_mat(I2,iN))),'r--','LineWidth',2');
    plot([sigma/2 sigma/2],[ 0 15],'--k')
    s_c = 2*acosh(1/(1-exp(-sigma/2)));
    s_d = -log(1/(1+2*(cosh(sigma)-1)/N));
    s_c = 1/N^2
    plot([s_c s_c],[ 0 15],'--k')
    plot([s_d s_d],[ 0 15],'--g')
    
end
ylabel('Re(\lambda)');
xlabel('s');
%%
figure;
axes('FontSize',24);
hold on
for iN=1:N
    
    I1=imag(epsilon_s_mat(:,iN))==0;
    I2=imag(epsilon_s_mat(:,iN))~=0;
    hold on
    plot(S(I1),(imag(epsilon_s_mat(I1,iN))),'k-','LineWidth',2');
    plot(S(I2),(imag(epsilon_s_mat(I2,iN))),'r--','LineWidth',2');
end
ylabel('Im(\lambda)');
xlabel('s')
% print(gcf, '-depsc2', '/Users/danielhurowitz/PROJ/NEG/Figs/tbWell_imag.eps')
%% Inverse localization length
figure;
axes('FontSize',24);
hold on;
plot(E_vec,kappa/N);

plot(E_vec(sign(expKappa)>0),kappa(sign(expKappa)>0)/N,'g.-',...
    E_vec,ones(1,length(E_vec))*log(2*cosh(N*s/2)-2)/N,'r',...
    E_vec,ones(1,length(E_vec))*s/2,'r');
a=get(gca);
%  plot(E_vec,4*sigma./(E_vec),'k')
axis([a.XLim a.YLim])

xlabel('\lambda');
ylabel('\kappa');
title(['s=',num2str(s),', \sigma=',num2str(sigma)]);