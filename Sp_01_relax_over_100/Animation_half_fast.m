close all
clear
clc

load('it.mat')
% it = it -270;

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
load('alpha_s.mat')
load('n_nodes_const_per_spine.mat')
load('n_el_sing.mat')


n_smooth = 128;
tvec = -1:2/n_smooth:1;
phi_l = (tvec-1).*tvec/2;
phi_c = (tvec+1).*(1-tvec);
phi_r = (tvec+1).*tvec/2;

vidObj = VideoWriter('released_half_fast');
vidObj.FrameRate = 12;
vidObj.Quality = 100;
open(vidObj);

Vels = figure;

load('time_vec.mat')
    


for step = 0:5:it
    load(['u_',num2str(step),'.mat'])
    load(['w_',num2str(step),'.mat'])
    load(['spine_lengths_',num2str(step),'.mat'])
    %locating nodes
    [Nodes_rz,~] ...
        = nodes_relocator_split_v02(spine_lengths, ...
                                    spine_feet_locator, ...
                                    n_spines, ...
                                    n_spines_incr,n_v,alpha_s, ...
                                    n_nodes_const_per_spine, ...
                                    n_el_sing);
    figure(Vels);
    %Plotting vecotr field                            
    quiver(Nodes_rz(:,1),Nodes_rz(:,2),u,w)
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
    set(gca,'xlim',[0 2],'ylim',[0 2])
    grid on
    set(gcf, 'Position', get(0, 'Screensize'));
    title(['t = ',sprintf('%8.3f',T_dim*time_vec(step+1)*1E9),' $n$s'],...
          'interpreter','latex','Fontsize',48) 
    xlabel('$x/R$','interpreter','latex','Fontsize',32)
    ylabel('$\frac{z}{R}\ \ $','interpreter','latex','Fontsize',48,'Rotation',0)
    drawnow
    hold off
    currFrame = getframe(gcf);
    writeVideo(vidObj,currFrame);
     
end
close(vidObj)