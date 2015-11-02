%% Plots electrostatic potential and field lines
%Input: epsilon_j_mat(iS,:)
%       s_vec(iS)
%       g        
iS=3
s=s_vec(iS);

allCharges = epsilon_j_mat(iS,:);           %set charge locations
% maxCharge = 4;                             %truncate charge vector
maxCharge = 200;                             %truncate charge vector
charges=allCharges(allCharges<=maxCharge);  
delta_X=diff(charges);
charges=charges(1:end-1);

x_vec=[linspace(0,maxCharge,1000)];        %x values to evaluate V,E
y_vec=[linspace(0,maxCharge/2,1000)];       %y values to evaluate V,E

% y_vec=[linspace(-DX/2,-1e-2,1000),linspace(1e-2,DX/2,1000)];
[x,y]=meshgrid(x_vec,y_vec);
phi=-sum(log(g))/N;
E_x=0;
E_y=0;

%Calculate electrostatic potential and field lines
for iE=1:length(allCharges)
    X = allCharges(iE);
    phi = phi +  log((x-X).^2+y.^2)/2/N;
    E_x = E_x -(x-X)./((x-X).^2+y.^2)/N;
    E_y = E_y - y./((x-X).^2+y.^2)/N;
end
%
%     phi_x=0;
%     E_x_x=0;
%     for iE=1:length(charges)
%         X = charges(iE);
%         phi_x = phi_x +  log((x_vec-X).^2)/2/N;
%         E_x_x = E_x_x -(x_vec-X)./((x_vec-X).^2)/N;
% %         E_y = E_y - y./((x-X).^2)/N;
%     end


%%
%Set coordinates to begin drawing field lines from
startx = charges;%+delta_X/2;
startx =linspace(charges(1),maxCharge,10);
starty1 = zeros(size(startx))+1e-15;
starty2 = zeros(size(startx))-1e-15;

figure;
axes('FontSize',24);
hold on;
imagesc(x_vec,y_vec,phi,[s/2-.3, s/2+.3]);colormap(flipud(RdBu)/255)
    imagesc(x_vec,-y_vec,(phi))
%     imagesc(x_vec,y_vec,phi,[1.4 3.2]);colormap(jet)
%     imagesc(x_vec,-y_vec,(phi),[1.4 3.2])
%     plot(charges,zeros(size(charges)),'g*')
%     plot(charges,zeros(size(charges)),'k.')
h=plot(charges(charges<=maxCharge),zeros(size(charges(charges<=maxCharge))),'k.','MarkerSize',40)
set(h,'Color',[0.3 0.3 0.3])
% set(h,'Color',[.5 0 .5])
% set(gch,'Color',[0.1 0.1 0.1])
% plot(startx(1),0,'*r')
% [CC,h_c1]=contour(x,y,phi,10,'k') ;
% [CC,h_c2]=contour(x,-y,phi,10,'k') ;
% set(h_c1,'Color',[1 0.2 0])
% set(h_c2,'Color',[1 0.2 0])

s1=streamline(x,y,-E_x,-E_y,startx,starty1);
set(s1,'Color','k');

s2=streamline(x,-y,-E_x,E_y,startx,starty2);
set(s2,'Color','k');

% 
% s1=streamline(x,y,-E_x,-E_y,charges(1),1e-15);
% set(s1,'Color','green');
% s2=streamline(x,y,-E_x,-E_y,charges(1),-1e-15);
% 
% s2=streamline(x,-y,-E_x,E_y,charges(1),-1e-15);
% set(s2,'Color','green');



% [C,h]=contour(x,y,phi,[s/2 s/2],'r','LineWidth',15);
% set(h,'Color',[1 1 1])
% contour(x,y,phi,[log(2*(cosh(s/2)-1)) log(2*(cosh(s/2)-1))],'r','LineWidth',4);
% [C,h]=contour(x,-y,phi,[s/2 s/2],'r','LineWidth',15);
% set(h,'Color',[1 1 1])
h=plot(-epsilon_s_mat(iS,real(-epsilon_s_mat(iS,:))<maxCharge),'.','MarkerSize',20)
set(h,'Color',[0 .7 0])
% plot(startx,zeros(size(startx)),'c*')
% axis([0 4.2 -2.1 2.1])
% axis square
xlabel('Re[\lambda]')
ylabel('Im[\lambda]')
axis([0 maxCharge -maxCharge/2 maxCharge/2 ])
% axis([0 32 -16 16 ])
% axis square
% title(['s=',num2str(s)])
%     hold off
colorbar
% colorbar
%     print(gcf, '-depsc2', ['/Users/danielhurowitz/PROJ/NEG/Figs/electrostatics_',num2str(iS),'.eps'])
%%
% x=linspace(0,100,1e3);
%     pause
%     pause(1)
% end
% [C,h]= contour(x,y,phi,[s/2 s/2]);
x_c=fliplr(C(1,2:end));
y_c=fliplr(C(2,2:end));
x_i=linspace(x_c(1,1),x_c(1,end),1000);
y_i=interp1(x_c,y_c,x_i);
% x_i=0;y_i=0;
E_x_0=0;
E_y_0=0;
phi_0=0;
for iE=1:length(allCharges)
    X = allCharges(iE);
    phi_0 = phi_0 +  log((x_i-X).^2+y_i.^2)/2/N;
    E_x_0 = E_x_0 -(x_i-X)./((x_i-X).^2+y_i.^2)/N;
    E_y_0 = E_y_0 - y_i./((x_i-X).^2+y_i.^2)/N;
end

plot(x_i,sqrt(x_i*4*exp(s/2)*cosh(sigma/2)/sinh(sigma)),'g','LineWidth',2)
plot(x_i,-sqrt(x_i*4*exp(s/2)*cosh(sigma/2)/sinh(sigma)),'g','LineWidth',2)


delta = sqrt(diff(x_i).^2+diff(y_i).^2);
N_lambda = cumsum(sqrt(E_x(1:end-1).^2+E_y(1:end-1).^2).*delta);
figure;plot(x_i,[0,N_lambda],sort(allCharges),(1:N)/N)