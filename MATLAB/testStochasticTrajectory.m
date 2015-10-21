% Generates stochastic trajectory and plots it on a ring

%% Set input parameters
clear all;
N=5000                  %Number of sites

randInd1 = randperm(N); %randomizes field
randInd2 = randperm(N); %randomizes couplings
alpha=Inf;0.7;0.1;      %P(w) ~ w^(alpha-1)
sigma=5;                %width of stochastic field distribution

%% Initialize model
s_1_2=log(sinh(sigma/2)/(sigma/2))*2;
s_1=log(sinh(sigma)/sigma);
s_2=log(sinh(2*sigma)/(2*sigma))/2;

x_n = (1:N)/N;
x_n = x_n(randInd2);
g=(x_n./N).^(1/alpha);

x0 = linspace(-1,1,N);
x0 = x0(randInd1) - mean(x0);
%% 
colors='rgbkm'
s_vec=[s_1_2*0.5,s_1*1.1]; %affinity values 
t_end=5000;         %End time of simulation
r0=(round(rand(1)*N)); %Initial position
for iS=1:length(s_vec)
    s=s_vec(iS)
    x=x0*sigma+s;
    
    
    w_p = [g(1:N-1).*exp(x(1:N-1)/2),g(N)*exp(-x(N)/2)];
    w_m = [g(1:N-1).*exp(-x(1:N-1)/2),g(N)*exp(x(N)/2)];
    %     for iR=1:1
    
    [x_vec,t_vec] = generateStochasticTrajectory(w_m,w_p,r0,t_end);
    %         x_mat{iR}=x_vec;
    %         t_mat{iR}=t_vec;
    %     end
    %
    % figure;
    % hold on;
    % cellfun(@plot,t_mat,x_mat);
    [P(iS,:),Pa(iS,:),I(iS),er]=ness(w_p,w_m,0);
    figure;
    
    polar((mod(x_vec,N)/N)*2*pi,t_vec,'r')%,colors(iS))
    hold on
    
    h= polar(2*pi*(-1:N-2)/N,t_end*(1.1+3*P(iS,:)),'b');
%     title(['s=',num2str(s)])
    set(h,'LineWidth',4);
    % grid on
    %     end
    set(gca,'FontSize',24)
    axis image
end

%%
% W = generateRateMatrix(x,g);
% p0=zeros(1,N);
% p0=ones(1,N)/N;
% % p0(1)=1;
% t_vec2 =linspace(0,1e5,1e2);
% figure;
% 
% 
% for iT=1:length(t_vec2)
%     polar(2*pi*(-1:N-2)/N,1+P,'k');
%     hold on
%     t=t_vec2(iT);
%     p(:,iT) = expm(W*t)*p0';
%     polar(2*pi*(0:N-1)/N,1+p(:,iT)')
%     %     plot(1:N,p(:,iT),'.')
%     %     axis([1 N 0 0.1])
%     pause;
%     hold off
% end
% %%
% for iT=1:3000
%     plot(1:N,p(:,iT),'.')
% end
