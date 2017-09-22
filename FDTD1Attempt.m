clear all;
close all;

eps0= 8.854E-12; mu0 = 4*pi*1e-7; c = 3e8; 
fsource = 3e8; 
xmax = 5; xmin = 0; dx = 0.01;
dt = dx / c;
x = (xmin:dx:xmax);

sigma = 0.0*(heaviside(x-2.5));
mu = mu0*ones(1,length(x));
eps = eps0*ones(1,length(x));
rho = zeros(1,length(x));

H_z = zeros(1,length(x));
E_y = zeros(1,length(x));

C_a = 2*eps - dt*sigma;
C_b = 2*dt/dx;
C_denom = 1./(2*eps + dt*sigma);

D_a = 2*mu-rho*dt;
D_b = 2*dt/(dx);
D_denom = 1./(2*mu + rho*dt);
z = zeros(length(x));
play = true;
n = 1;
while play
   E_y(250) = sin(-2*pi*fsource*dt*n);%E_y(250)+2*(1-exp(-((n-1)/50)^2))*sin(-2*pi*fsource*dt*n);
%    if n == 1
%        E_y(250) = 0.5;
%    end
   
   H_z(1) = H_z(2); 
   
   % This is H_z(i,n+1)
   for i=2:length(x)-1 
      H_z(i) = D_a(i)*D_denom(i)*H_z(i)-D_b*D_denom(i)*(E_y(i+1)-E_y(i)); 
   end
   
   E_y(length(x)) = E_y(length(x)-1); 
   
   for i=2:length(x)-1 
      E_y(i) = C_a(i)*C_denom(i)*E_y(i)-C_b*C_denom(i)*(H_z(i)-H_z(i-1));
   end
   
   figure(1)
   hold off
   plot(x,E_y,'b');
%    plot3(x,E_y,z,x,z,377*H_z)
%    axis([3*dx xmax -2 2 -2 2]);
   grid on;
   title('E (blue) and 377*H (red)');
   hold on
   plot(x,377*H_z,'r');
   xlabel('x (m)');
   ylabel('y ');
   axis([3*dx xmax -2 2])
%    zlabel('z (H,m)');
   pause(1e-3);
   n = n+1
end



