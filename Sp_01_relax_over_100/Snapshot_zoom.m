close all
clear
clc

width = .005;
height = .5E-2;
scale = .01;

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
load('n_spines_incr.mat')
load('alpha_s.mat')
load('n_nodes_const_per_spine.mat')
load('n_el_sing.mat')
load('U_cl_vec.mat')


n_smooth = 128;
tvec = -1:2/n_smooth:1;
phi_l = (tvec-1).*tvec/2;
phi_c = (tvec+1).*(1-tvec);
phi_r = (tvec+1).*tvec/2;

% vidObj = VideoWriter('released_half');
% vidObj.FrameRate = 5;
% vidObj.Quality = 100;
% open(vidObj);

Vels = figure;

load('delta_t.mat')
% load('it_start_hs1.mat')
% tvec = delta_t*(0:it_start_hs-1);
% tvec = delta_t*(0:it);
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
load('time_vec.mat')
    
labels = ['a','b','c','d'];
count = 0 ;
for step = floor([it/20 it/2 3*it/4 it-2])
    count = count + 1;
    load(['u_',num2str(step),'.mat'])
    load(['w_',num2str(step),'.mat'])
    load(['spine_lengths_',num2str(step),'.mat'])
    %locating nodes
    [Nodes_rz,~] ...
        = nodes_relocator_split_v02(spine_lengths, ...
                                     spine_feet_locator,n_spines, ...
                                     n_spines_incr,n_v,alpha_s, ...
                                     n_nodes_const_per_spine, ...
                                     n_el_sing);
    figure(Vels);
    h_x = Nodes_rz(1,1);
    h_y = Nodes_rz(1,2);
    %Plotting vecotr field                            
    quiver(Nodes_rz(:,1)-h_x,Nodes_rz(:,2)-h_y,u-U_cl_vec(step),w,scale)
    hold on  
    %Plotting boundary 1
    for ii = 1:n1_el
        xvec = Nodes_rz(l_1(ii,1),1)*phi_l ...
               + Nodes_rz(l_1(ii,2),1)*phi_c ...
               + Nodes_rz(l_1(ii,3),1)*phi_r;
        yvec = Nodes_rz(l_1(ii,1),2)*phi_l ...
               + Nodes_rz(l_1(ii,2),2)*phi_c ...
               + Nodes_rz(l_1(ii,3),2)*phi_r;
        plot(xvec-h_x,yvec-h_y,'b','LineWidth',4)
    end
    %Plotting boundary 2
    for ii = 1:n2_el
        xvec = Nodes_rz(l_2(ii,1),1)*phi_l ...
               + Nodes_rz(l_2(ii,2),1)*phi_c ...
               + Nodes_rz(l_2(ii,3),1)*phi_r;
        yvec = Nodes_rz(l_2(ii,1),2)*phi_l ...
               + Nodes_rz(l_2(ii,2),2)*phi_c ...
               + Nodes_rz(l_2(ii,3),2)*phi_r;
        plot(xvec-h_x,yvec-h_y,'g','LineWidth',4)
    end
    %Plotting boundary 3
    for ii = 1:n3_el
        xvec = Nodes_rz(l_3(ii,1),1)*phi_l ...
               + Nodes_rz(l_3(ii,2),1)*phi_c ...
               + Nodes_rz(l_3(ii,3),1)*phi_r;
        yvec = Nodes_rz(l_3(ii,1),2)*phi_l ...
               + Nodes_rz(l_3(ii,2),2)*phi_c ...
               + Nodes_rz(l_3(ii,3),2)*phi_r;
        plot(xvec-h_x,yvec-h_y,'r','LineWidth',4)
    end
    set(gca,'FontSize',24)
    axis equal
    set(gca,'xlim',[-.75*width .25*width],'ylim',[0 height],'FontSize',24)%,'xtick',0:.2:2)
    text(-.25,2.05,[labels(count),')'],'FontSize',24)
    grid on
    set(gcf, 'Position', get(0, 'Screensize'));
    title(['t = ',sprintf('%8.3f',T_dim*time_vec(step+1)*1E9),' $n$s'],...
          'interpreter','latex','Fontsize',48) 
    xlabel('$x/R$','interpreter','latex','Fontsize',32)
    ylabel('$\frac{z}{R}\ \ $','interpreter','latex','Fontsize',48,'Rotation',0)
    drawnow
    hold off
%     currFrame = getframe(gcf);
%     writeVideo(vidObj,currFrame);
    print('-depsc',['Spreading_zoom',num2str(count),'.eps'])
end
% close(vidObj)