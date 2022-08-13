close all
clear
clc

load('it.mat')
load('spine_feet_locator.mat')
load('n_spines.mat')
load('n_spines_incr.mat')
load('n_v_start_far.mat')
load('n_v.mat')
load('n_v_near.mat')
load('n1_el_near.mat')
load('n1_el.mat')
load('l_1.mat')
load('n2_el_near.mat')
load('n2_el.mat')
load('l_2.mat')
load('n3_el.mat')
load('l_3.mat')
load('n4_el.mat')
load('l_4.mat')
load('U_cl.mat')


% load('spine_lengths_ini.mat')
% spine_lengths = spine_lengths_ini;
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


Vels = figure;
 
bound_to_glob1_map = zeros(1,n_spines);
bound_to_glob1_map(1) = 1;
kk = 1;
for ee = 1:n1_el
    kk = kk+1;
    bound_to_glob1_map(kk) = l_1(ee,2);
    kk = kk+1;
    bound_to_glob1_map(kk) = l_1(ee,3);
end

bound_to_glob2_map = zeros(1,n_spines);
bound_to_glob2_map(1) = 1;
kk = 1;
for ee = 1:n1_el
    kk = kk+1;
    bound_to_glob2_map(kk) = l_2(ee,2);
    kk = kk+1;
    bound_to_glob2_map(kk) = l_2(ee,3);
end

for step = it
    load(['us1_',num2str(step),'.mat'])
    load(['ws1_',num2str(step),'.mat'])
    load(['us2_',num2str(step),'.mat'])
    load(['ws2_',num2str(step),'.mat'])
    load(['spine_lengths_',num2str(step),'.mat'])
    %locating nodes
    [Nodes_rz,~] ...
        = nodes_relocator_split_v01(spine_lengths, ...
                                    spine_feet_locator, ...
                                    n_spines, ...
                                    n_spines_incr,n_v);
    figure(Vels);
    %Plotting vecotr field 
    quiver(Nodes_rz(bound_to_glob1_map,1), ...
           Nodes_rz(bound_to_glob1_map,2), ...
           us1-U_cl(step+1), ...
           ws1, ...
           .01)
    hold on  
    quiver(Nodes_rz(bound_to_glob2_map,1), ...
           Nodes_rz(bound_to_glob2_map,2), ...
           us2-U_cl(step+1), ...
           ws2, ...
           .1)
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
%     %Plotting boundary 4
%     for ii = 1:n4_el
%         xvec = Nodes_rz(l_4(ii,1),1)*phi_l ...
%                + Nodes_rz(l_4(ii,2),1)*phi_c ...
%                + Nodes_rz(l_4(ii,3),1)*phi_r;
%         yvec = Nodes_rz(l_4(ii,1),2)*phi_l ...
%                + Nodes_rz(l_4(ii,2),2)*phi_c ...
%                + Nodes_rz(l_4(ii,3),2)*phi_r;
%         plot(xvec,yvec,'r','LineWidth',2)
%     end
    set(gca,'FontSize',16)
    axis equal
    set(gca,'xlim',[-.001 .005]+Nodes_rz(1,1),'ylim',[-0.001 .002])
    grid on
    set(gcf, 'Position', get(0, 'Screensize'));
    title(['time step = ',num2str(step)]) 
    drawnow
    hold off  
    pause
end