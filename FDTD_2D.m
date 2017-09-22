clear all;
close all;

eps0= 8.854E-12; mu0 = 4*pi*1e-7; c = 3e8; 
fsource = 3e8; 
xmax = 5; xmin = 0; dx = 0.01;
dt = dx / c;
x = (xmin:dx:xmax);
y = (xmin:dx:xmax);

sigma = zeros(length(y),length(x));
mu = mu0*ones(length(y),length(x));
eps = eps0*ones(length(y),length(x));
rho = zeros(length(y),length(x));

H_z = zeros(length(y),length(x));
E_y = zeros(length(y),length(x));
E_x = zeros(length(y),length(x));

C_a = 2*eps - dt*sigma;
C_b = 2*dt/dx;
C_denom = 1./(2*eps + dt*sigma);

D_a = 2*mu-rho*dt;
D_b = 2*dt/(dx);
D_denom = 1./(2*mu + rho*dt);
z = zeros(length(x));
play = true;
n = 1;
while n < 60
   E_y(250,250) = 2*(1-exp(-((n-1)/50)^2))*sin(-2*pi*fsource*dt*n);%E_y(250,250)+
   H_z(:,1) = H_z(:,2); 
   H_z(1,:) = H_z(2,:);
   
   % This is H_z(i,n+1)
   for i=2:length(x)-1
       for j = 2:length(y) -1
           H_z(j,i) = D_a(j,i).*D_denom(j,i).*H_z(j,i)-D_b*D_denom(j,i).*(E_y(j,i+1)-E_y(j,i))+D_b*D_denom(j,i).*(E_x(j,i)-E_x(j+1,i)); 
       end
   end
   
   E_y(:,length(x)) = E_y(:,length(x)-1); 
   E_y(length(y),:) = E_y(length(y)-1,:);
   E_x(:,length(x)) = E_x(:,length(x)-1); 
   
   E_x(length(y),:) = E_x(length(y)-1,:);
   
   for i=2:length(x)-1
       for j = 2:length(y)-1
           E_y(j,i) = C_a(j,i).*C_denom(j,i).*E_y(j,i)-C_b*C_denom(j,i).*(H_z(j,i)-H_z(j,i-1));
           E_x(j,i) = C_a(j,i).*C_denom(j,i).*E_x(j,i)-C_b*C_denom(j,i).*(H_z(j+1,i)-H_z(j,i));
       end
   end
   
   figure(1)
   hold off
%    surf(x,y,E_y)
   plot(x,E_y(250,:),'b');
%    plot3(x,E_y(4,:),z,x,z,377*H_z(4,:))
%    axis([3*dx xmax -2 2 -2 2]);
   grid on;
   title('E (blue) and 377*H (red)');
   hold on
   plot(x,377*H_z(250,:),'r');
   xlabel('x (m)');
   ylabel('y ');
%    axis([3*dx xmax -2 2])
%    colorbar
%    zlabel('z (H,m)');
   pause(1e0);
   n = n+1
end
