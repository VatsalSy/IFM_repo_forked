close all
clear
clc

width = .1;
height = 5E-2;
scale = .3;

load('it.mat')
% it = it - 270;


load('spine_feet_locator.mat')
load('n_spines.mat')
load('n_spines_incr.mat')
load('n_v.mat')
load('n1_el.mat')
load('l_1.mat')
load('n2_el.mat')
load('l_2.mat')
load('n3_el.mat')
load('l_3.mat')
load('T_dim.mat')
load('U_cl_vec.mat')
load('hs.mat')
load('tvec.mat')

load(['u_',num2str(0),'.mat'],'u')
load(['w_',num2str(0),'.mat'],'w')

load('spine_lengths_ini.mat')
spine_lengths = spine_lengths_ini;
save('spine_lengths0.mat','spine_lengths')
% [Nodes_rz,~] ...
%         = nodes_relocator_split_v00(spine_lengths, ...
%                                     spine_feet_locator, ...
%                                     n_spines, ...
%                                     n_spines_incr,n_v);
%for plotting
n_smooth = 128;
t_vec = -1:2/n_smooth:1;
phi_l = (t_vec-1).*t_vec/2;
phi_c = (t_vec+1).*(1-t_vec);
phi_r = (t_vec+1).*t_vec/2;

vidObj = VideoWriter('released_half_zoom_moving_frame');
vidObj.FrameRate = 5;
vidObj.Quality = 100;
open(vidObj);

Vels = figure;

 


for step = 0:1:it
    load(['u_',num2str(step),'.mat'])
    load(['w_',num2str(step),'.mat'])
    load(['spine_lengths_',num2str(step),'.mat'])
    %locating nodes
    [Nodes_rz,~] ...
        = nodes_relocator_split_v01(spine_lengths, ...
                                    spine_feet_locator, ...
                                    n_spines, ...
                                    n_spines_incr,n_v);
    figure(Vels);
    %Plotting vecotr field                            
    quiver(Nodes_rz(:,1),Nodes_rz(:,2),u-U_cl_vec(step+1),w,scale)
    hold on  
    %Plotting boundary 1
    for ii = 1:n1_el
        xvec = Nodes_rz(l_1(ii,1),1)*phi_l ...
               + Nodes_rz(l_1(ii,2),1)*phi_c ...
               + Nodes_rz(l_1(ii,3),1)*phi_r;
        yvec = Nodes_rz(l_1(ii,1),2)*phi_l ...
               + Nodes_rz(l_1(ii,2),2)*phi_c ...
               + Nodes_rz(l_1(ii,3),2)*phi_r;
        plot(xvec,yvec,'b','LineWidth',2)
    end
    %Plotting boundary 2
    for ii = 1:n2_el
        xvec = Nodes_rz(l_2(ii,1),1)*phi_l ...
               + Nodes_rz(l_2(ii,2),1)*phi_c ...
               + Nodes_rz(l_2(ii,3),1)*phi_r;
        yvec = Nodes_rz(l_2(ii,1),2)*phi_l ...
               + Nodes_rz(l_2(ii,2),2)*phi_c ...
               + Nodes_rz(l_2(ii,3),2)*phi_r;
        plot(xvec,yvec,'g','LineWidth',2)
    end
    %Plotting boundary 3
    for ii = 1:n3_el
        xvec = Nodes_rz(l_3(ii,1),1)*phi_l ...
               + Nodes_rz(l_3(ii,2),1)*phi_c ...
               + Nodes_rz(l_3(ii,3),1)*phi_r;
        yvec = Nodes_rz(l_3(ii,1),2)*phi_l ...
               + Nodes_rz(l_3(ii,2),2)*phi_c ...
               + Nodes_rz(l_3(ii,3),2)*phi_r;
        plot(xvec,yvec,'r','LineWidth',2)
    end
    set(gca,'FontSize',16)
    axis equal
    set(gca,'xlim',[Nodes_rz(1,1)-width/2 Nodes_rz(1,1)+width/2],'ylim',[0 height])
    grid on
    set(gcf, 'Position', get(0, 'Screensize'));
    title(['t = ',sprintf('%9.4f',T_dim*tvec(step+1)*1E9),' $n$s'],...
          'interpreter','latex','Fontsize',48) 
    xlabel('$x/R$','interpreter','latex','Fontsize',32)
    ylabel('$\frac{z}{R}\ \ $','interpreter','latex','Fontsize',48,'Rotation',0) 
    drawnow
    hold off
    currFrame = getframe(gcf);
    writeVideo(vidObj,currFrame);
%     pause
end
close(vidObj)