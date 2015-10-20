clear all

N=11;
n=0:N-1;


s_vec=linspace(0,5,1000);
x0=zeros(1,N);
g=ones(1,N);
% x0=ones(1,N);
for iS = 1:length(s_vec)
    
    s = s_vec(iS);
    x = x0 + s;
%     g=exp(B0);
%     g=ones(1,N);
    lambda_s(iS,:) = SolveNHR(x,g);
    lambda(iS,:) =2*( cosh(s/2)*(cos(2*pi*n/N)-1)+1i*sinh(s/2)*sin(2*pi*n/N));
    S(iS)=sum(x)/N;
end

figure;
axes('FontSize',24);
hold on
plot(lambda_s,'b','LineWidth',4);
plot(lambda,'r--','LineWidth',1);
