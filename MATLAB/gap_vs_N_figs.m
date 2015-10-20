figure;
axes('FontSize',24);
hold on

load gap_500
N=500;
gap=epsilon_s_mat(:,N-1);
I=imag(gap)~=0;
plot(s_vec(imag(gap)==0),N^2*abs(real(gap(imag(gap)==0))),'g','LineWidth',4);

load gap_1000
N=1000;
gap=epsilon_s_mat(:,N-1);
I=imag(gap)~=0;
plot(s_vec(imag(gap)==0),N^2*abs(real(gap(imag(gap)==0))),'r','LineWidth',4);

load('disorder_N_3000_full.mat')
gap=epsilon_s_mat(:,N-1);
% gap2=epsilon_s_mat(:,N-2);
I=imag(gap)~=0;

plot(s_vec(imag(gap)==0),N^2*abs(real(gap(imag(gap)==0))),'LineWidth',4);


load gap_500
N=500;
gap=epsilon_s_mat(:,N-1);
I=imag(gap)~=0;
I=I(10:end);
plot(s_vec(I),N^2*abs(real(gap(I))),'g--','LineWidth',4);


load gap_1000
N=1000;
gap=epsilon_s_mat(:,N-1);
I=imag(gap)~=0;
I=I(4:end);
plot(s_vec(I),N^2*abs(real(gap(I))),'r--','LineWidth',4);


load('disorder_N_3000_full.mat')
gap=epsilon_s_mat(:,N-1);
% gap2=epsilon_s_mat(:,N-2);
I=imag(gap)~=0;
I=I(4:end);
plot(s_vec(I),N^2*abs(real(gap(I))),'--','LineWidth',4);

plot(s_vec,4*pi^2*cosh(s_vec/2),'k','LineWidth',2);
plot(s_vec,2*pi^2*exp(s_vec/2-s_1_2)*cosh(sigma/2),'c','LineWidth',2)

plot(s_vec([10,20,27,40]),-N^2*real(gap([10,20,27,40])),'gd','MarkerFaceColor','g','LineWidth',2)


plot([s_1_2 s_1_2],[0 800],'--k','LineWidth',2)
plot([s_1 s_1],[0 800],'--k','LineWidth',2)
plot([s_2 s_2],[0 800],'--k','LineWidth',2)
plot([sigma sigma],[0 800],'--k','LineWidth',2)
xlabel('s');
ylabel('N^2\Delta_{\lambda}');
legend(['N=500 ';...
        'N=1000';...
        'N=3000'])
axis([0 7 0 800]);
% print(gcf, '-depsc2', ['/Users/danielhurowitz/PROJ/NEG/Figs/gap_s_N.eps'])
%%
figure;
axes('FontSize',24);
hold on

load gap_500
N=500;
gap=epsilon_s_mat(:,N-1);
I=imag(gap)~=0;
plot(s_vec(imag(gap)==0),N*abs(real(gap(imag(gap)==0))),'g','LineWidth',4);

load gap_1000
N=1000;
gap=epsilon_s_mat(:,N-1);
I=imag(gap)~=0;
plot(s_vec(imag(gap)==0),N*abs(real(gap(imag(gap)==0))),'r','LineWidth',4);

load('disorder_N_3000_full.mat')
gap=epsilon_s_mat(:,N-1);
% gap2=epsilon_s_mat(:,N-2);
I=imag(gap)~=0;

plot(s_vec(imag(gap)==0),N*abs(real(gap(imag(gap)==0))),'LineWidth',4);


load gap_500
N=500;
gap=epsilon_s_mat(:,N-1);
I=imag(gap)~=0;
I=I(10:end);
plot(s_vec(I),N*abs(real(gap(I))),'g--','LineWidth',4);


load gap_1000
N=1000;
gap=epsilon_s_mat(:,N-1);
I=imag(gap)~=0;
I=I(4:end);
plot(s_vec(I),N*abs(real(gap(I))),'r--','LineWidth',4);


load('disorder_N_3000_full.mat')
gap=epsilon_s_mat(:,N-1);
% gap2=epsilon_s_mat(:,N-2);
I=imag(gap)~=0;
I=I(4:end);
plot(s_vec(I),N*abs(real(gap(I))),'--','LineWidth',4);

% plot(s_vec,4*pi^2*cosh(s_vec/2),'k','LineWidth',2);

plot(s_vec([10,20,27,40]),-N*real(gap([10,20,27,40])),'gd','MarkerFaceColor','g','LineWidth',2)


plot([s_1_2 s_1_2],[0 800],'--k','LineWidth',2)
plot([s_1 s_1],[0 800],'--k','LineWidth',2)
plot([s_2 s_2],[0 800],'--k','LineWidth',2)
plot([sigma sigma],[0 800],'--k','LineWidth',2)
xlabel('s');
ylabel('N\Delta_{\lambda}');
legend(['N=500 ';...
        'N=1000';...
        'N=3000'])
axis([0 7 0 1]);
% print(gcf, '-depsc2', ['/Users/danielhurowitz/PROJ/NEG/Figs/gap_s_N.eps'])
%%

figure;
axes('FontSize',24);
hold on
grid on
% plot(s_vec(imag(gap)==0),abs(real(gap(imag(gap)==0))),'b','LineWidth',4);
% plot(s_vec(imag(gap)~=0),abs(real(gap(imag(gap)~=0))),'--b','LineWidth',4);
plot(s_vec,Diffusion*4*pi^2*N^2,'g','LineWidth',4);
plot(s_vec,abs(real(gap))*N^2,'--','LineWidth',4);
% plot(s_vec,1./(rho_0*N),'--r')
plot(s_vec,4*pi^2*cosh(s_vec/2),'k','LineWidth',2);

% plot(s_vec([10,20,27,40]),-N*real(gap([10,20,27,40])),'gd','MarkerFaceColor','g','LineWidth',2)
% plot(s_vec,1/sigma*(((s_vec+sigma)/2).^2-((s_vec-sigma)/2).^2)/2/N,'g')
% plot(s_vec,1./(exp(-s_vec/2).*sinh(sigma/2)/sigma),'g');
% plot(s_vec,sigma^2*pi^2*exp(s_vec/2)/sinh(sigma/2)/2/N^2,'r');
% plot(s_vec,sigma^2*pi^2*exp(s_vec/2)*coth(sigma/2)/sinh(sigma/2)/N^2/2,'g')
plot(s_vec,2*pi^2*exp(s_vec/2-s_1_2)*cosh(sigma/2),'r','LineWidth',2)
% plot(s_vec,2*pi^2*exp((s_vec/2-s_1_2))*tanh(sigma/2)/N^2,'c')

plot([s_1_2 s_1_2],[0 N^2*max(abs(real(gap)))],'--k','LineWidth',2)
plot([s_1 s_1],[0 N^2*max(abs(real(gap)))],'--k','LineWidth',2)
plot([s_2 s_2],[0 N^2*max(abs(real(gap)))],'--k','LineWidth',2)
plot([sigma sigma],[0 N^2*max(abs(real(gap)))],'--k','LineWidth',2)
xlabel('s');
ylabel('N^2\Gamma');
legend({'4\pi^2D';...
        'N^2 Re[\lambda_1]';...
        'Clean ring';...
        'Log-box'})
axis([0 7 0 700]);
% axtype(3)