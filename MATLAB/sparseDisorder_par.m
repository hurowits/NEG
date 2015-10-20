function  sparseDisorder_par(processNumber)
mkdir(['NEG/datasets/prcs',num2str(processNumber)]);
% addpath(genpath('/PROJ/'))
%% Initialization


N=1000;
sigma=5;
M_vec = [0:100:N-10,N-9:N-1];
% M_vec=498
% M_vec=N-2


randInd1 = randperm(N);
randInd2 = randperm(N);
% sigma=5;%sigma_vec(iSigma);;
g=ones(1,N);

s_1_2=log(sinh(sigma/2)/(sigma/2))*2;
s_1=log(sinh(sigma)/sigma);
s_2=log(sinh(2*sigma)/(2*sigma))/2;


%             s_vec=2*sigma;
s_vec = linspace(0,s_1_2,500);
numComplex = zeros(length(M_vec),length(s_vec));

x0 = linspace(-1,1,N);

x0 = x0(randInd1) - mean(x0);

for iM=1:length(M_vec)
    
    M=M_vec(iM);
    randInd3 = randInd2(1:M);
    x0_1 = x0;
    x0_1(randInd3)=0;
    %     if (max(randInd3)==N)
    %         [c,I]=max(randInd3);
    %         randInd3(I)=1;
    %     end
    %     x0_1(randInd3+1)=0;
    
    for iS=1:length(s_vec)
        %                 mu=mu_vec(iS);
        
        s = s_vec(iS);
        %     s=4;
        x=x0_1*sigma+s;
        s_numeric = mean(x);
        
        [epsilon_s] = SolveNHR(x,g);  %find true eigenvalues
        epsilon_s_mat(iS,:)=epsilon_s;
        numComplex(iM,iS)=sum(imag(epsilon_s)~=0); %how many complex eigenvalues?
        
        gap(iS)=epsilon_s(N-1);
        
        if (numComplex(iM,iS)>10)
            s_c(iM)=mean(x);
            break;
        end
        
    end
end
    save (['NEG/datasets/prcs',num2str(processNumber),'/NEG_Dataset',num2str(iR),'.mat']);

% save(['s_c_',num2str(iR),'.mat']);
