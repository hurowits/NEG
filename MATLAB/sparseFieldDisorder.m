%% Initialization
clear all
iP=6;

colors='rgbckymrgbckym'

options = optimset('TolX',1e-11,'TolFun',1e-11);

N=100;
% 
% M_vec = [0:50:N-50,N-9:N-1];
% M_vec = [0:10:N-10,N-9:N-1];
M_vec = [1:9,round(logspace(1,2,11))];[1:9,10:5:N];
% M_vec=round(logspace(log10(0.002),0,19)*N);
% M_vec=498
% M_vec=N-2
for iR = 1:100
    clc
    display(['realization number = ',num2str(iR)])
    tic ;
    randInd1 = randperm(N);
    randInd2 = randperm(N);
    sigma=5;%sigma_vec(iSigma);;
    g=ones(1,N);
    
    s_1_2=log(sinh(sigma/2)/(sigma/2))*2;
    s_1=log(sinh(sigma)/sigma);
    s_2=log(sinh(2*sigma)/(2*sigma))/2;
    
    
    %             s_vec=2*sigma;
    s_vec = linspace(0,s_1_2,500);
    numComplex = zeros(length(M_vec),length(s_vec));
    
    x0 = linspace(-1,1,N);
    
    x0 = x0(randInd1) - mean(x0);
    %     hwb=waitbar(0,'Please wait');
    x0_1=zeros(1,N);
    for iM=1:length(M_vec)
        %         waitbar(iM/length(M_vec))
%         iM;
        M=M_vec(iM);
        randInd3 = randInd2(1:M);
        
        x0_1(randInd3)= x0(randInd3);
        
        %         x0_1 = x0;
        %         x0_1(randInd3)=0;
        %         %     if (max(randInd3)==N)
        %         [c,I]=max(randInd3);
        %         randInd3(I)=1;
        %     end
        %     x0_1(randInd3+1)=0;
        
        for iS=1:length(s_vec)
            %                 mu=mu_vec(iS);
            
            s = s_vec(iS);
            %     s=4;
            x=x0_1*sigma+s;
            %             s_numeric = mean(x);
            
            [epsilon_s] = SolveNHR(x,g);  %find true eigenvalues
            epsilon_s_mat(iS,:)=epsilon_s;
            numComplex(iM,iS)=sum(imag(epsilon_s)~=0); %how many complex eigenvalues?
            
            gap(iS)=epsilon_s(N-1);
            
            if (numComplex(iM,iS)>10)
                s_c(iM)=mean(x);
                s_c_mat(iR,iM)=s_c(iM);
                break;
            end
            
        end
    end
    toc
    %     save(['s_c_',num2str(iP),'_',num2str(iR),'.mat']);
    
end

% save(['s_c_',num2str(iP),'_',num2str(iR),'.mat']);

% close(hwb)
%%
%%
% iP=1;
% iN=1;
% for iP=4:4
%     for iR=1:10
%         % if (iN==11)continue;end
%         load(['s_c_',num2str(iP),'_',num2str(iR),'.mat']);
%         s_c_mat2(iN,:)=s_c;
%         iN=iN+1;
%     end
% end
% s_c_mat=[s_c_mat(1:10,:);s_c_mat(12:20,:)]
%%
figure;
axes('FontSize',24);
hold on;
grid on

max(abs(s_c_mat),1)

plot(M_vec/N,(abs(s_c_mat))/s_1_2,'b.','MarkerSize',10)
y=mean(abs(s_c_mat),1)/s_1_2;
y(end)=1;
plot(M_vec/N,y,'rs','MarkerSize',10,'MarkerFaceColor','r')

% plot((N-M_vec)/N,sigma/N/s_1_2*(N-M_vec),'--k','LineWidth',4)
plot(M_vec/N,sqrt(M_vec)*sigma/N/s_1_2/2,'-k','LineWidth',2)
% axtype(3);
% 
% upper = mean(abs(s_c_mat),1)/s_1_2 + std(abs(s_c_mat/s_1_2),1)/2;
% lower = mean(abs(s_c_mat),1)/s_1_2 - std(abs(s_c_mat/s_1_2),1)/2;
axtype(3)
% hf=ciplot(lower,upper,M_vec/N,'r')

% set(hf,'FaceAlpha',0.5)
% errorbar((N-M_vec)/N,abs(s_c)/s_1_2,diff(s_vec(1:2))*ones(1,length(M_vec)),'-*','LineWidth',4)
% errorbar((N-M_vec)/N,mean(abs(s_c_mat),1)/s_1_2,std(abs(s_c_mat),1)/s_1_2,'b.','LineWidth',4)
% plot((N-M_vec)/N,(abs(s_c_mat))/s_1_2,'bo','LineWidth',2)
% plot((N-M_vec)/N,mean(abs(s_c_mat),1)/s_1_2,'r*','LineWidth',4)
% 
% 
% % plot((N-M_vec)/N,sigma/N/s_1_2*(N-M_vec),'--k','LineWidth',4)
% plot((N-M_vec)/N,sqrt(N-M_vec)*sigma/N/s_1_2/2,'-k','LineWidth',2)
xlabel('$M/N$','Interpreter','Latex');
ylabel('$s_c/s_{1/2}$','Interpreter','Latex')
% % axis([0 1 0 1])
% print(gcf, '-depsc2', ['/Users/danielhurowitz/PROJ/NEG/Figs/s_c_sparse_100_loglog.eps'])

% save('s_c_sparseFieldDisorder');