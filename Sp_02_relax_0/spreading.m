close all
clear
clc

dtheta = pi/1000;

figure
theta_vec = 0:dtheta:2*pi;
x_circle = cos(theta_vec);
z_circle = sin(theta_vec);
plot(x_circle,z_circle,'k','LineWidth',2)
hold on
plot([-1 1],[-1.1 -1.1],'k','LineWidth',2)
quiver(0,.5,0,-1,'k','LineWidth',2,'MaxHeadSize',.4)
axis equal
axis off
print('-depsc','figure1.eps')
pause

figure
angle = .9*pi/2;
theta_vec = -angle:dtheta:pi+angle;
x_circle = cos(theta_vec);
z_circle = sin(theta_vec);
plot(x_circle,z_circle,'k','LineWidth',2)
hold on
plot([-1 1],[-sin(angle) -sin(angle)],'k','LineWidth',2)
quiver([0 .5 -.5],[.75 -.5 -.5],[0 .15 -.15],[-.5 -.15 -.15],'k','LineWidth',2,'MaxHeadSize',.8)
axis equal
axis off
print('-depsc','figure2.eps')
pause

figure
angle = .7*pi/2;
alpha = 50;
beta = .1;
theta_vec = -angle:dtheta:pi+angle;
x_circle = cos(theta_vec);
z_circle = sin(theta_vec);
plot(x_circle,z_circle,'k','LineWidth',2)
hold on
plot([-1 1],[-sin(angle) -sin(angle)],'k','LineWidth',2)
quiver([.5 -.5],[-.75 -.75],[.25 -.25],[0 0],'k','LineWidth',2,'MaxHeadSize',.8)
x_solid = 0:.01:cos(angle)/2;
x_solid = [fliplr(-x_solid),x_solid(2:end)];
z_solid = beta*exp(-alpha*x_solid.^2)-sin(angle);
fill(x_solid,z_solid,'r')
axis equal
axis off
print('-depsc','figure3.eps')
pause

figure
angle = .7*pi/2;
alpha = 15;
beta = .1;
theta_vec = -angle:dtheta:pi+angle;
x_circle = cos(theta_vec);
z_circle = sin(theta_vec);
plot(x_circle,z_circle,'k','LineWidth',2)
hold on
plot([-1 1],[-sin(angle) -sin(angle)],'k','LineWidth',2)
x_solid = 0:.01:cos(angle);
x_solid = [fliplr(-x_solid),x_solid(2:end)];
z_solid = beta*exp(-alpha*x_solid.^2)-sin(angle);
fill(x_solid,z_solid,'r')
axis equal
axis off
print('-depsc','figure4.eps')
pause

figure
angle = .7*pi/2;
theta_vec = -angle:dtheta:pi+angle;
x_circle = cos(theta_vec);
z_circle = sin(theta_vec);
plot(x_circle,z_circle,'k','LineWidth',2)
hold on
plot([-1 1],[-sin(angle) -sin(angle)],'k','LineWidth',2)
x_solid = [0 cos(angle) cos(angle-dtheta:-dtheta:pi/4) cos(pi/4:dtheta:pi/2)];
x_solid = [fliplr(-x_solid),x_solid(2:end)];
z_solid = [-sin(angle) -sin(angle) -sin(angle-dtheta:-dtheta:pi/4) sin(pi/4:dtheta:pi/2)-2*sin(pi/4)];
z_solid = [fliplr(z_solid),z_solid(2:end)];
fill(x_solid,z_solid,'r')
axis equal
axis off
print('-depsc','figure5.eps')
pause

figure
angle = .7*pi/2;
theta_vec = -angle:dtheta:pi+angle;
x_circle = cos(theta_vec);
z_circle = sin(theta_vec);
plot(x_circle,z_circle,'k','LineWidth',2)
hold on
plot([-1 1],[-sin(angle) -sin(angle)],'k','LineWidth',2)
x_solid = [0 cos(angle) cos(angle-dtheta:-dtheta:-pi/2)];
x_solid = [fliplr(-x_solid),x_solid(2:end)];
z_solid = [-sin(angle) -sin(angle) -sin(angle-dtheta:-dtheta:-pi/2)];
z_solid = [fliplr(z_solid),z_solid(2:end)];
fill(x_solid,z_solid,'r')
axis equal
axis off
print('-depsc','figure6.eps')