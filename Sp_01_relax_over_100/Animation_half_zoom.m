close all
clear
clc

width = .1;
height = 5E-2;
scale = .1;

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
load('tvec.mat')


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
tvec = -1:2/n_smooth:1;
phi_l = (tvec-1).*tvec/2;
phi_c = (tvec+1).*(1-tvec);
phi_r = (tvec+1).*tvec/2;

vidObj = VideoWriter('released_half_zoom');
vidObj.FrameRate = 5;
vidObj.Quality = 100;
open(vidObj);

Vels = figure;
load('delta_t.mat')
% load('it_start_hs1.mat')
load('it.mat')
% tvec = delta_t*(0:it_start_hs-1);
tvec = delta_t*(0:it);
% count = it_start_hs;
% for cc = 2:9
%     load(['delta_t_hs_',num2str(cc-1),'.mat'])
%     load(['it_start_hs',num2str(cc),'.mat'])
%     tvec = [tvec,tvec(end)+delta_t*(1:(it_start_hs-count))];
%     count = it_start_hs;
% end
% load(['delta_t_hs_',num2str(cc),'.mat'])
% % load('it.mat')
% tvec = [tvec,tvec(end)+delta_t*(1:it)];
 


for step = 0:10:it
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
    quiver(Nodes_rz(:,1),Nodes_rz(:,2),u,w,scale)
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
    title(['t = ',sprintf('%8.5f',T_dim*tvec(step+1)*1E9),' $n$s'],...
          'interpreter','latex','Fontsize',48) 
    xlabel('$x/R$','interpreter','latex','Fontsize',32)
    ylabel('$\frac{z}{R}\ \ $','interpreter','latex','Fontsize',48,'Rotation',0) 
    drawnow
    hold off
    currFrame = getframe(gcf);
    writeVideo(vidObj,currFrame);
     
end
close(vidObj)