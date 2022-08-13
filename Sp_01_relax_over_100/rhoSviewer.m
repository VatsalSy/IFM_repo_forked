close all
clear 
clc

load('it.mat')
% it = it-270;
load('n_spines.mat')
load('Ds.mat')
load('Dg.mat')
load('n1_el.mat')
load('l_2.mat')
load('l_1.mat')
load('spine_feet_locator.mat')
load('n_spines_incr.mat')
load('n_v.mat')
load('T_dim.mat')
load('time_vec.mat')
load('U_cl_vec.mat')
load('sigma_1e_dim.mat')
load('sigma_2e_dim.mat')
load('alpha_1.mat')
load('theta_c_vec.mat')
load('alpha_s.mat')
load('n_nodes_const_per_spine.mat')
load('n_el_sing.mat')

rhos2_mat = zeros(n_spines,it+1);
sigma2_mat = zeros(n_spines,it+1);
for ii = 0:it
    load(['rhos2_',num2str(ii),'.mat'])
    rhos2_mat(:,ii+1) = rhos2;
    load(['sigma2_',num2str(ii),'.mat'])
    sigma2_mat(:,ii+1) = sigma2;
end
maximum = max(max(rhos2_mat));
minimum = min(min(rhos2_mat));
range = maximum-minimum;


bound_to_glob2_map = zeros(1,n_spines);
bound_to_glob2_map(1) = 1;
kk = 1;
for ee = 1:n1_el
    kk = kk+1;
    bound_to_glob2_map(kk) = l_2(ee,2);
    kk = kk+1;
    bound_to_glob2_map(kk) = l_2(ee,3);
end
save('bound_to_glob2_map.mat','bound_to_glob2_map')

bound_to_glob1_map = zeros(1,n_spines);
bound_to_glob1_map(1) = 1;
kk = 1;
for ee = 1:n1_el
    kk = kk+1;
    bound_to_glob1_map(kk) = l_1(ee,2);
    kk = kk+1;
    bound_to_glob1_map(kk) = l_1(ee,3);
end
save('bound_to_glob1_map.mat','bound_to_glob1_map')

load('hs.mat')    
%Derivatives with respect to xi, used to find the tangent at the joint
phi1_xi = zeros(3,2);
xi(1) = -1;
xi(2) = 0;
xi(3) = 1;
for pp = 1:3
    phi1_xi(1,pp) = xi(pp)-.5;
    phi1_xi(2,pp) = -2*xi(pp);
    phi1_xi(3,pp) = xi(pp)+.5;
end
%Mesh velocity on boundary 1
mesh_r_vel_b1 = zeros(n_spines,it+1);
mesh_z_vel_b1 = zeros(n_spines,it+1);
%Mesh velocity on boundary 2
mesh_r_vel_b2 = zeros(n_spines,it+1);
mesh_z_vel_b2 = zeros(n_spines,it+1);
load(['spine_lengths_',num2str(0),'.mat'])
[Nodes_rz_minus1,~] ...
    = nodes_relocator_split_v02(spine_lengths, ...
                                spine_feet_locator, ...
                                n_spines, ...
                                n_spines_incr,n_v,alpha_s, ...
                                n_nodes_const_per_spine, ...
                                n_el_sing);

    
load(['spine_lengths_',num2str(ii),'.mat'])
[Nodes_rz,~] ...
    = nodes_relocator_split_v02(spine_lengths, ...
                            spine_feet_locator, ...
                            n_spines, ...
                            n_spines_incr,n_v,alpha_s, ...
                            n_nodes_const_per_spine, ...
                            n_el_sing);
delta_t = time_vec(2)-time_vec(1);
mesh_r_vel_b1(:,2) = (  Nodes_rz(bound_to_glob1_map,1) ...
                         - Nodes_rz_minus1(bound_to_glob1_map,1) ...
                        )/delta_t;
mesh_z_vel_b1(:,2) = (  Nodes_rz(bound_to_glob1_map,2) ...
                         - Nodes_rz_minus1(bound_to_glob1_map,2) ...
                        )/delta_t;
mesh_r_vel_b2(:,2) = (  Nodes_rz(bound_to_glob2_map,1) ...
                         - Nodes_rz_minus1(bound_to_glob2_map,1) ...
                        )/delta_t;
mesh_z_vel_b2(:,2) = (  Nodes_rz(bound_to_glob2_map,2) ...
                         - Nodes_rz_minus1(bound_to_glob2_map,2) ...
                        )/delta_t;
for ii = 2:it
    delta_t_prev = delta_t;
    delta_t = time_vec(ii+1)-time_vec(ii);
    Nodes_rz_minus2 = Nodes_rz_minus1;
    Nodes_rz_minus1 = Nodes_rz;
    load(['spine_lengths_',num2str(ii),'.mat'])
    [Nodes_rz,~] ...
        = nodes_relocator_split_v02(spine_lengths, ...
                                spine_feet_locator, ...
                                n_spines, ...
                                n_spines_incr,n_v,alpha_s, ...
                                n_nodes_const_per_spine, ...
                                n_el_sing);
    ratio_n = delta_t/delta_t_prev;
    a_n = (1+2*ratio_n)/(1+ratio_n);
    a_n_minus1 = -(1+ratio_n);
    a_n_minus2 = ratio_n^2/(1+ratio_n);
    mesh_r_vel_b1(:,ii+1) ...
            = (  a_n*Nodes_rz(bound_to_glob1_map,1) ...
               + a_n_minus1*Nodes_rz_minus1(bound_to_glob1_map,1) ...
               + a_n_minus2*Nodes_rz_minus2(bound_to_glob1_map,1) ...
              )/(delta_t);
        mesh_z_vel_b1(:,ii+1) ...
            = (  a_n*Nodes_rz(bound_to_glob1_map,2) ...
               + a_n_minus1*Nodes_rz_minus1(bound_to_glob1_map,2) ...
               + a_n_minus2*Nodes_rz_minus2(bound_to_glob1_map,2) ...
              )/(delta_t);
        mesh_r_vel_b2(:,ii+1) ...
            = (  a_n*Nodes_rz(bound_to_glob2_map,1) ...
               + a_n_minus1*Nodes_rz_minus1(bound_to_glob2_map,1) ...
               + a_n_minus2*Nodes_rz_minus2(bound_to_glob2_map,1) ...
              )/(delta_t);
        mesh_z_vel_b2(:,ii+1) ...
            = (  a_n*Nodes_rz(bound_to_glob2_map,2) ...
               + a_n_minus1*Nodes_rz_minus1(bound_to_glob2_map,2) ...
               + a_n_minus2*Nodes_rz_minus2(bound_to_glob2_map,2) ...
              )/(delta_t);   
end

    

flow_into_b2 = zeros(1,it+1);
Dens2 = figure('units','normalized','outerposition',[0 0 1 1]);

vidObj = VideoWriter('Surface2_variables');
vidObj.FrameRate = 5;
vidObj.Quality = 100;
open(vidObj);

for step = 0:1:it
    load(['spine_lengths_',num2str(step),'.mat'])
    [Nodes_rz,~] ...
        = nodes_relocator_split_v02(spine_lengths, ...
                                    spine_feet_locator, ...
                                    n_spines, ...
                                    n_spines_incr,n_v,alpha_s, ...
                                    n_nodes_const_per_spine, ...
                                    n_el_sing);
    plot(Nodes_rz(bound_to_glob2_map,1),rhos2_mat(:,step+1),'b', ...
         'LineWidth',2)
    hold on
    plot([min(Nodes_rz(:,1)) max(Nodes_rz(bound_to_glob2_map,1))], ...
         [Ds Ds],'--b','LineWidth',2)
    plot(Nodes_rz(bound_to_glob2_map,1),sigma2_mat(:,step+1),'k', ...
         'LineWidth',2)
    plot([min(Nodes_rz(:,1)) max(Nodes_rz(bound_to_glob2_map,1))], ...
         [sigma_2e_dim/sigma_1e_dim sigma_2e_dim/sigma_1e_dim],'--k','LineWidth',2)
    load(['us2_',num2str(step),'.mat'])    
    plot(Nodes_rz(bound_to_glob2_map,1),-(us2-mesh_r_vel_b2(:,step+1)),'-x','LineWidth',2)
    flow_into_b2(step+1) = -rhos2_mat(1,step+1)*(us2(1)-mesh_r_vel_b2(1,step+1));
    load(['w_',num2str(step),'.mat'])
    plot(Nodes_rz(bound_to_glob2_map,1),w(bound_to_glob2_map),'-+','LineWidth',2)
    set(gca,'FontSize',16)%,'ylim',[minimum-range/10 maximum+range/10]
    xlabel('$x/R$','interpreter','latex','FontSize',24)
%     ylabel('$\frac{\rho^s_2}{\rho^s_0}\ \ \ \ \ \ \ \ $', ...
%            'interpreter','latex','rotation',0,'FontSize',32)
    title(['t = ',sprintf('%9.2f',T_dim*time_vec(step+1)*1E9),' $n$s'],...
          'interpreter','latex','Fontsize',36) 
    grid on
    legend('$\rho^{s_2}/\rho^{s}_0$','$\rho^{s_2}_e/\rho^s_0$', ...
           '$\sigma^2/\sigma^1_e$','$\sigma^2_e/\sigma^1_e$', ...
           '$u^s_2-u_c$','$w$','interpreter','latex', ...
           'location','West','FontSize',24)
    drawnow
    pause(.05)
    hold off
    currFrame = getframe(gcf);
    writeVideo(vidObj,currFrame);
end
close(vidObj)

rhos1_mat = zeros(n_spines,it+1);
sigma1_mat = zeros(n_spines,it+1);
for ii = 0:it
    load(['rhos1_',num2str(ii),'.mat'])
    rhos1_mat(:,ii+1) = rhos1;
    load(['sigma1_',num2str(ii),'.mat'])
    sigma1_mat(:,ii+1) = sigma1;
end
maximum = max(max(max(rhos1_mat)),1E-9);
minimum = min(min(min(rhos1_mat)),-1E-9);
range = maximum-minimum;

%Surface1
Dens1 = figure('units','normalized','outerposition',[0 0 1 1]);
   
flow_into_b1 = zeros(1,it+1);
flow_into_b1_prime = zeros(1,it+1);
% tangent_b1 = figure;

vidObj = VideoWriter('Surface1_variables');
vidObj.FrameRate = 5;
vidObj.Quality = 100;
open(vidObj);

for step = 0:it
%     figure(Dens1)
    load(['spine_lengths_',num2str(step),'.mat'])
    [Nodes_rz,~] ...
        = nodes_relocator_split_v02(spine_lengths, ...
                                    spine_feet_locator, ...
                                    n_spines, ...
                                    n_spines_incr,n_v,alpha_s, ...
                                    n_nodes_const_per_spine, ...
                                    n_el_sing);
    plot(Nodes_rz(bound_to_glob1_map,1),rhos1_mat(:,step+1),'b', ...
         'LineWidth',2)
    hold on
    plot([min(Nodes_rz(:,1)) max(Nodes_rz(bound_to_glob1_map,1))], ...
         [Dg Dg],'--b','LineWidth',2)
    plot(Nodes_rz(bound_to_glob1_map,1),sigma1_mat(:,step+1),'k', ...
         'LineWidth',2)
    plot([min(Nodes_rz(:,1)) max(Nodes_rz(bound_to_glob1_map,1))], ...
         [1 1],'--k','LineWidth',2)
    %finding tangent vectors
    ee = 1;
    t1_r = zeros(n_spines,1);
    t1_z = zeros(n_spines,1);
    t1_r(1) = (phi1_xi(:,1)')*Nodes_rz(l_1(ee,:),1);
    t1_z(1) = (phi1_xi(:,1)')*Nodes_rz(l_1(ee,:),2);
    t1_r(2) = (phi1_xi(:,2)')*Nodes_rz(l_1(ee,:),1);
    t1_z(2) = (phi1_xi(:,2)')*Nodes_rz(l_1(ee,:),2);
    t1_r(3) = .5*(phi1_xi(:,3)')*Nodes_rz(l_1(ee,:),1);
    t1_z(3) = .5*(phi1_xi(:,3)')*Nodes_rz(l_1(ee,:),2);
    kk = 3;
    for ee = 2:n1_el
        t1_r(kk) = t1_r(kk) + .5*(phi1_xi(:,1)')*Nodes_rz(l_1(ee,:),1);
        t1_z(kk) = t1_z(kk) + .5*(phi1_xi(:,1)')*Nodes_rz(l_1(ee,:),2);
        kk = kk+1;
        t1_r(kk) = (phi1_xi(:,2)')*Nodes_rz(l_1(ee,:),1);
        t1_z(kk) = (phi1_xi(:,2)')*Nodes_rz(l_1(ee,:),2);
        kk = kk+1;
        t1_r(kk) = .5*(phi1_xi(:,3)')*Nodes_rz(l_1(ee,:),1);
        t1_z(kk) = .5*(phi1_xi(:,3)')*Nodes_rz(l_1(ee,:),2);
    end
    t1_r(kk) = t1_r(kk) + .5*(phi1_xi(:,3)')*Nodes_rz(l_1(ee,:),1);
    t1_z(kk) = t1_z(kk) + .5*(phi1_xi(:,3)')*Nodes_rz(l_1(ee,:),2);
    for ii = 1:n_spines
        factor = max(abs(t1_r(ii)),abs(t1_z(ii)));
        t1_r(ii) = t1_r(ii)/factor;
        t1_z(ii) = t1_z(ii)/factor;
    end
    angle = atan2(t1_z(1),t1_r(1));
    co = cos(angle);
    si = sin(angle);
    t1_r = t1_r./sqrt(t1_r.^2+t1_z.^2);
    t1_z = t1_z./sqrt(t1_r.^2+t1_z.^2);
    load(['us1_',num2str(step),'.mat'])
    load(['ws1_',num2str(step),'.mat'])
    flow_into_b1(step+1) = rhos1_mat(1,step+1) ...
                           *(- co ...%t1_r(1) ...
                               .*(us1(1)-mesh_r_vel_b1(1,step+1)) ...
                             + si ...%t1_z(1) ...
                               .*(ws1(1)-mesh_z_vel_b1(1,step+1)));
    flow_into_b1_prime(step+1) = rhos1_mat(1,step+1) ...
                                 *(- cos(theta_c_vec(step+1)) ...
                                     .*(us1(1)-mesh_r_vel_b1(1,step+1)) ...
                                   + sin(theta_c_vec(step+1)) ...
                                     .*(ws1(1)-mesh_z_vel_b1(1,step+1)) ...
                                  );
    set(gca,'FontSize',16)%,'ylim',[minimum-range/10 maximum+range/10]
    xlabel('$x/R$','interpreter','latex','FontSize',24)
%     ylabel('$\frac{\rho^s_1}{\rho^s_0}\ \ \ \ \ \ \ \ $', ...
%            'interpreter','latex','rotation',0,'FontSize',32)
    title(['t = ',sprintf('%9.2f',T_dim*time_vec(step+1)*1E9),' $n$s'],...
          'interpreter','latex','Fontsize',36) 
    grid on
    
    plot(Nodes_rz(bound_to_glob1_map,1), ...
         (  t1_r.*(us1-mesh_r_vel_b1(:,step+1)) ...
          + t1_z.*(ws1-mesh_z_vel_b1(:,step+1)))/10,'-x')
    load(['u_',num2str(step),'.mat'])
    load(['w_',num2str(step),'.mat'])
    plot(Nodes_rz(bound_to_glob1_map,1), ...
         (- alpha_1*t1_z.*(u(bound_to_glob1_map)-us1) ...
          + alpha_1*t1_r.*(w(bound_to_glob1_map)-ws1))/10,'-+')
    legend('$\rho^{s_1}/\rho^s_0$','$\rho^{s_1}_e/\rho^s_0$', ...
           '$\sigma^1/\sigma^1_e$','$\sigma^1_e/\sigma^1_e$', ...
           '$(u^s_{1\|}-c_{\|})$','$u_n-u^s_n$', ...
           'interpreter','latex','location','West', ...
           'FontSize',24)
    drawnow
    hold off
    
%     figure(tangent_b1)
%     quiver(Nodes_rz(bound_to_glob1_map,1), ...
%            Nodes_rz(bound_to_glob1_map,2), ...
%            t1_r, ...
%            t1_z)
%     axis equal
    pause(.05)
    currFrame = getframe(gcf);
    writeVideo(vidObj,currFrame);
end
close(vidObj)

figure
plot(time_vec(1:it+1)*T_dim*1E9,-flow_into_b1_prime,'LineWidth',6)
hold on
plot(time_vec(1:it+1)*T_dim*1E9,flow_into_b2,'color',[0 1 0],'LineWidth',2)
% plot(time_vec(1:it+1)*T_dim*1E9,-flow_into_b1,'--','LineWidth',6)
xlabel('T [ns]','FontSize',24)
set(gca,'FontSize',24)
grid on
legend('Mass flux out of free surf','Mass flux into liquid-solid surf', ...
       'location','SouthEast','FontSize',24)