clear all

N=100;


% 
% x0 = 2*linspace(0,1,N) - 1;
% B0 = 2*linspace(0,1,N) - 1;
% randInd1 = randperm(N);
% randInd2 = randperm(N);
% 
% x0 = x0(randInd1) - mean(x0);
% B0 = B0(randInd2) - mean(B0);
% 
% 
% x0 = zeros(1,N);zeros(1,N);
% n = floor(N/2);
% s1=10;
% x0(n)=x0(n)-s1;
% x0(n+1)=x0(n+1)+s1;
% 
% x0(7)=0;
% x0(8)=0;

%     x0 = (-1).^(0:N-1);
s_vec=linspace(0,20,10000);
x0=zeros(1,N);
g=ones(1,N);
g_0=1;1e-8;.0001;
g(1)=g_0;
sigma=1;2*pi;5*pi;5*pi;
% g(1)=0;1e-8;.0001;1e-8;
% x0=[-ones(1,N/2),ones(1,N/2)];
for iS = 1:length(s_vec)
    
    s = s_vec(iS);
    x = x0*sigma + s;
%     g=exp(B0);
%     g=ones(1,N);
    lambda_s(iS,:) = SolveNHR(x,g);
    
    S(iS)=sum(x)/N;
end

gap = -lambda_s(:,end-1);

I = abs(imag(gap))<1e-5;
s_c=log(sinh(sigma*2)/sigma/2);
figure;plot(S,log10(real(lambda_s)),'.')
figure;plot(S,(imag(lambda_s)),'-')
figure;plot(lambda_s,'.')
%%




figure(1);
% axes('FontSize',24)
hold on;
plot(s_vec(I),real(gap(I)),'r','LineWidth',4);
plot(s_vec(~I),real(gap(~I)),'r--','LineWidth',4);
plot(s_vec,4*pi^2/N^2*ones(size(s_vec)),'--k')
plot(s_vec,s_vec.^2/4,'--k')
plot(s_vec,(s_vec-sigma).^2/4,'--k')

% plot(s_vec,-real(lambda_s),'--');

% plot([s_c s_c],[0 0.03],'--k','LineWidth',2)
xlabel('s');
ylabel('Gap');
% title(['\sigma=',num2str(sigma),', g=',num2str(g_0)]);
% 
% legend(['I - real     ';...
%         'I - complex  ';...
%         'II - real    ';...
%         'II - complex ';...
%         'III - real   ';...
%         'III - complex']);
%     axis([0 8 -0.1 15])
%%
I=imag(gap)==0;
figure(3);
% axes('FontSize',24)
hold on;
% grid on
plot(s_vec,-real(lambda_s),'b-','LineWidth',1);
plot(s_vec,s_vec.^2/4,'k')
plot(s_vec(I),(sigma-s_vec(I)).^2/4,'--r','LineWidth',2);
plot(s_vec(~I),4*pi^2/N^2*ones(1,length(s_vec(~I))) ,'.r','LineWidth',2)
% plot(s_vec,4*s_vec/N+pi^2/N^2,'k--','LineWidth',4)

% axis([2 7 -0.05 0.5])
xlabel('s');
ylabel('Re[\lambda]')
%%
f=cos(0.5*sqrt(x.^2-sigma^2/2))-(x.^2+sigma^2/4).*sin(0.5*sqrt(x.^2-sigma^2/2))./sqrt(x.^2-sigma^2/2);