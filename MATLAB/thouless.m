%% Thouless formula for Anderson model
N = 5;
sigma =0;.1;3;
s=1;
gamma=2*ones(1,N);((rand(1,N)-0.5)*sigma); 
% gamma=exp(s+(rand(1,N)-0.5)*sigma);
% gamma = gamma-mean(gamma);
% gamma = gamma-mean(gamma);
W = diag(gamma,0) +...
    diag(-ones(1,N-1),1) +...
    diag(-ones(1,N-1),-1);


W(1,N)=-1;0.001;
W(N,1)=-1;0.001;
[v,lambda_j] = eig(W);
lambda_j=diag(lambda_j)
% figure;hist(lambda_j,100)

E_vec = linspace(min(lambda_j),max(lambda_j),1e4);
for iE = 1:length(E_vec)
    E=E_vec(iE);
%     kappa(iE) = prod(E-lambda_j);
     kappa(iE)=sum(log(abs(E-lambda_j)))/N;
end



figure;
axes('FontSize',24);
hold on
% plot(E_vec,log(abs(kappa)));
plot(-E_vec,((kappa)),'LineWidth',2);
a=get(gca);
% plot(E_vec,ones(1,length(E_vec))*log(2*cosh(s*N/2)-2)/N);
% plot(E_vec,sigma^2/24./(4-(E_vec-2*sinh(sigma/2)/sigma).^2),'--r','LineWidth',2)
plot(-E_vec,sigma^2/24./(4-(E_vec-mean(gamma)).^2),'--r','LineWidth',2)

axis([a.XLim a.YLim]);
xlabel('E');
ylabel('\kappa = \xi^{-1}');
legend(['Thouless';...
       'FGR     '])
% axis([-3 3 -0.3 0.3])
% print(gcf, '-depsc2', '/Users/danielhurowitz/PROJ/NEG/Figs/anderson_box2.eps')