sigma=5
E_vec=linspace(0,100,1e3);
mu_vec=linspace(0.001,3,3);
[mu,E]=meshgrid(mu_vec,E_vec);
s = (1./mu).*log(sinh(sigma*mu)./(sigma*mu));
% s=sigma*mu/2;
D0 =sqrt(2*(cosh(s)-1)./s.^2);
x = sqrt((2304*D0.^3/sigma^4).*E);
N = sigma^2./(24*pi^2*D0.^2.*(besselj(mu,x).^2+bessely(mu,x).^2));

figure;
contour(mu,E,N,[1 1],'k','LineWidth',2)
xlabel('\mu');
ylabel('E');
title('N(E)');

figure;
contour(s,E,N,[1 1])
xlabel('s');
ylabel('E');
title('N(E)');

