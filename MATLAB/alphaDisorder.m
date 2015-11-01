%Find s_c vs. alpga
%% Initialization
clear all
iP=6;
alpha_vec=0.1:0.1:1.5;
sigma=0
colors='rgbckymrgbckym'

options = optimset('TolX',1e-11,'TolFun',1e-11);

N=500;
%
M_vec=[10,50,100,500,1000];
for iM=1:length(M_vec)
    N=M_vec(iM);
    for iA=1:length(alpha_vec)
        alpha = alpha_vec(iA)
        for iR = 1:20
            
            display(['M=',num2str(N),' alpha=',num2str(alpha),' realization number = ',num2str(iR)])
            tic ;
            randInd1 = randperm(N);
            randInd2 = randperm(N);
            
            g=ones(1,N);
            
            x_n = (1:N)/N;
            x_n = x_n(randInd2);
            g=(x_n./N).^(1/alpha);
            
            s_vec = linspace(0,1,200);
            numComplex = zeros(1,length(s_vec));
            
            x0 = linspace(-1,1,N);
            
            x0 = x0(randInd1) - mean(x0);
            %     hwb=waitbar(0,'Please wait');
            
            
            
            
            for iS=1:length(s_vec)
                %                 mu=mu_vec(iS);
                
                s = s_vec(iS);
                %     s=4;
                x=x0*sigma+s;
                %             s_numeric = mean(x);
                
                [epsilon_s] = SolveNHR(x,g);  %find true eigenvalues
%                 epsilon_s_mat(iS,:)=epsilon_s;
                numComplex(iM,iA,iS)=sum(imag(epsilon_s)~=0); %how many complex eigenvalues?
                
                gap(iS)=epsilon_s(N-1);
                
                if (numComplex(iM,iA,iS)>1)
                    s_c=mean(x);
                    s_c_mat(iM,iR,iA)=s_c;
                    break;
                end
                
            end
        end
        toc
        %     save(['s_c_',num2str(iP),'_',num2str(iR),'.mat']);
        
    end
end
% save(['s_c_',num2str(iP),'_',num2str(iR),'.mat']);

% close(hwb)
%%
%%
% iP=1;
% iN=1;
% for iP=1:3
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

plot(alpha_vec,(abs(s_c_mat)),'b.','MarkerSize',10)
y=mean(abs(s_c_mat),1);
% y(end)=1;
plot(alpha_vec,y,'rs','MarkerSize',10,'MarkerFaceColor','r')

% plot((N-M_vec)/N,sigma/N/s_1_2*(N-M_vec),'--k','LineWidth',4)
% plot(M_vec/N,sqrt(M_vec)*sigma/N/s_1_2/2,'-k','LineWidth',2)
% axtype(3);
%
% upper = mean(abs(s_c_mat),1)/s_1_2 + std(abs(s_c_mat/s_1_2),1)/2;
% lower = mean(abs(s_c_mat),1)/s_1_2 - std(abs(s_c_mat/s_1_2),1)/2;
% axtype(3)
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
xlabel('\alpha');
ylabel('s_c')
% % axis([0 1 0 1])
% print(gcf, '-depsc2', ['/Users/danielhurowitz/PROJ/NEG/Figs/s_c_sparse_100_loglog.eps'])