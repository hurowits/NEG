%Find s_c vs. alpga
%% Initialization
% clear all
iP=6;
alpha_vec=[0.1:0.1:3];
sigma=0

options = optimset('TolX',1e-11,'TolFun',1e-11);
numRealizations=100;
N=500;
%
M_vec=[10,100,500];%,500];
% s_vec = linspace(0,1,200);
% s_vec=logspace(-5,1,200);
% save('data/InitializationData_s_c_Vs_alpha')
for iM=1:length(M_vec)
    
    N=M_vec(iM);
    s_vec=[linspace(0,2,1e2)];
    
    epsilon_s_mat=zeros(length(s_vec),N);
    for iA=15%16:length(alpha_vec)
        alpha = alpha_vec(iA);
        for iR = 1:numRealizations
            
            display(['M=',num2str(N),' alpha=',num2str(alpha),' realization number = ',num2str(iR)])
            
            %             tic ;
            randInd1 = randperm(N);
            randInd2 = randperm(N);
            
            g=ones(1,N);
            
            x_n = (1:N)/N;
            x_n = x_n(randInd2);
            g=(x_n./N).^(1/alpha);
            
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
                epsilon_s_mat(iS,:)=epsilon_s;
                %                 epsilon_s_mat_ALL{iM,iA,iR,iS}=epsilon_s;
%                 numComplex(iM,iA,iR,iS)=sum(imag(epsilon_s)~=0); %how many complex eigenvalues?
                numComplex=sum(imag(epsilon_s)~=0); %how many complex eigenvalues?
                %                 save('data/numComplex');
                %                 gap(iS)=epsilon_s(N-1);
                
                if (numComplex>1)
                    s_c=mean(x);
                    s_c_mat(iM,iR,iA)=s_c;
                    break;
                end
                %
                %                                 elseif (numComplex(iM,iA,iS)>5)
                %                                     s_c5=mean(x);
                %                                     s_c5_mat(iM,iR,iA)=s_c5;
                %                                 elseif (numComplex(iM,iA,iS)>10)
                %                                     s_c10=mean(x);
                %                                     s_c10_mat(iM,iR,iA)=s_c10;
                %                                     break
                %                                 end
            end
%             save('data/alphaData3')
%             display(['saving epsilon_s_mat_N_',num2str(N),'_iAlpha_',num2str(iA),'_iR_',num2str(iR)])
            %             save(['data/epsilon_s_mat_N_',num2str(N),'_iAlpha_',num2str(iA),'_iR_',num2str(iR)],'epsilon_s_mat','randInd1','randInd2')
%             display('done')
        end
        %         toc;
        
        
        %     save(['s_c_',num2str(iP),'_',num2str(iR),'.mat']);
        
    end
    
end
% save(['s_c_',num2str(iP),'_',num2str(iR),'.mat']);

% close(hwb)
%%
%   numComplex(iM,iA,iS,iR)=sum(imag(epsilon_s)~=0); %how many complex eigenvalues?
% s_c_mat_MA=zeros(length(M_vec),length(alpha_vec));
% for iM=1:length(M_vec)
%     numComplex_M = squeeze(numComplex(iM,:,:,:));
%     for iA = 1:length(alpha_vec)
%         for iR = 1:numRealizations
%             for iS=1:length(s_vec)
%                 if (numComplex(iM,iA,iR,iS)>0)
%                     s_c_mat_FULL(iM,iA,iR)=s_vec(iS);
%                     iA
%                     iR
%                     iS
%                     break;
%                 end
%             end
%         end
%     end
% end

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
% s_c_mat2=s_c_mat;
% s_c_mat2(s_c_mat2==0)=Inf;
colors='gbr'

markerStyle='osd'
figure(1);
axes('FontSize',24);
hold on;
grid on
shiftFactor=0.02;
for iM=1:length(M_vec)
    M=M_vec(iM);
    s_c=squeeze(s_c_mat(iM,:,1:end))*M;
    %     max(abs(s_c_mat),1)
    y=mean(abs(s_c),1);
s_c2=s_c;        
s_c2(s_c==0)=Inf;
    y = sum(s_c2,1)./sum(s_c~=0,1);

    for iAlpha=1:length(alpha_vec)
        s_c_iA=s_c(:,iAlpha);
        y(iAlpha) = mean(s_c_iA(s_c_iA~=0));
    end
    
%     y(y==0)=NaN;
    plot(alpha_vec(1:end)+iM*shiftFactor,(abs(s_c2)),[colors(iM),'.'],'MarkerSize',10)
%     h(iM)=plot(alpha_vec(1:end)+iM*shiftFactor,y,[colors(iM),markerStyle(iM)],'MarkerSize',10,'MarkerFaceColor',colors(iM));
    h(iM)=plot(alpha_vec(1:end)+iM*shiftFactor,y,[colors(iM)],'LineWidth',4);
    % y(end)=1;
    
    
    
   
    xlabel('$\alpha$','Interpreter','Latex');
    ylabel('$Ns_c$','Interpreter','Latex')
    % % axis([0 1 0 1])
    % print(gcf, '-depsc2', ['/Users/danielhurowitz/PROJ/NEG/Figs/s_c_sparse_100_loglog.eps'])
end
legend(h(1:end),{['N=',num2str(M_vec(1))];['N=',num2str(M_vec(2))];['N=',num2str(M_vec(3))]},'Interpreter','Latex');
% plot(alpha_vec,exp(-2./(1+alpha_vec)))

% legend(h(1:end),{['N=',num2str(M_vec(2))];['N=',num2str(M_vec(3))]},'Interpreter','Latex');
axtype(2)
axis([0 1.5 0 100])
%%
[h1,x]=hist(log10(squeeze(s_c_mat(2:3,:,4))'.*repmat(M_vec(2:3),[numRealizations,1])))
handle=plot((-h1)/max(h1(:)),10.^x,'LineWidth',2)

%%
% [h1,x]=hist(log10(squeeze(s_c_mat2(3,:,3:6)).*M_vec(3)),10)
[h1,x]=hist(log10(squeeze(s_c_mat(3,:,3:6)).*M_vec(3)),10)
handle=plot((-h1)/max(h1(:)),10.^x,'--r','LineWidth',2)




set(handle(1),'Marker','.')
set(handle(2),'Marker','s')
set(handle(3),'Marker','d')
set(handle(4),'Marker','o')
% 
ah=axes('position',get(gca,'position'),...
            'visible','off','FontSize',20);
hleg=legend(ah,handle,[repmat('$\alpha =$',[4,1]),num2str(alpha_vec(3:6)')],'Location','SouthWest')
set(hleg,'Interpreter','Latex')

% [h,x]=hist((squeeze(s_c_mat2(1:3,:,3))'.*repmat(M_vec(1:3),[numRealizations,1])),10);
%%
[h,x]=hist(log10(squeeze(s_c_mat2(3,:,3:6)).*M_vec(3)),10)
figure;
axes('FontSize',24);
hold on;
plot((-h),x,'-o','LineWidth',4)
legend(num2str(alpha_vec(3:6)'))
ylabel('$\log(Ns_c)$','Interpreter','Latex');

% 
% figure(2);
% axes('FontSize',24);
% hold on;
% grid on
% [x,h]=hist(squeeze(M_vec(3)*s_c_mat(3,:,1:4)))