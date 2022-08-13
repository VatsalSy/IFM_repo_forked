clear
close all
clc

load('theta_c_ini.mat')
load('hs.mat','hs')
load('T_dim.mat')
load('Ds.mat')
load('Dg.mat')
load('it.mat')
load('alpha_s.mat')
load('n_nodes_const_per_spine.mat')
load('n_el_sing.mat')
% it = it-4;

load('time_vec.mat')
load('theta_c_vec.mat')
load('theta_m_vec.mat')
load('U_cl_vec.mat')


load('spine_feet_locator.mat')
load('n_spines.mat')
load('n_spines_incr.mat')
load('n_v.mat')



load('U_dim.mat')
load('theta_c_eq.mat')
rhos2_cl = zeros(1,it+1);
rhos1_cl = zeros(1,it+1);
for ii = 1:it+1
    load(['rhos2_',num2str(ii-1),'.mat'])
    rhos2_cl(ii) = rhos2(1);
    load(['rhos1_',num2str(ii-1),'.mat'])
    rhos1_cl(ii) = rhos1(1);    
end
load('theta_c_vec.mat','theta_c_vec')
figure
yyaxis left
theta_of_U = plot(U_cl_vec(1:it+1)*U_dim,theta_c_vec(1:it+1)*180/pi,'-b','linewidth',4)
hold on
theta_eq = plot([1.1*min(U_cl_vec)*U_dim 1.1*max(U_cl_vec)*U_dim],[theta_c_eq theta_c_eq]*180/pi,'--b','Linewidth',2)
set(gca,'xlim',[1.1*min(U_cl_vec(1:it+1))*U_dim 1.1*max(U_cl_vec(1:it+1))*U_dim], ...
        'ylim',180*[min(theta_c_eq,min(theta_c_vec(1:it+1)))-10/180 ...
                    max([theta_c_eq,theta_c_vec(1:it+1)])+10/180]/pi,'FontSize',32)
ylabel('$\theta_c [^o]\ \ \ \ $','interpreter','LaTeX','rotation',0,'FontSize',48)
xlabel('$U_{cl}[m/s]$','interpreter','latex','FontSize',32)
% legend([],'$\theta_c$','$\theta_{c,eq}$', ...
%        'interpreter','Latex','FontSize',16,'Location','SouthEast')
yyaxis right
rhos1_of_U = plot(U_cl_vec(1:it+1)*U_dim,rhos1_cl(1:it+1),'-k','LineWidth',2)
load('rhos_0_dim.mat')
load('rhos1_e_dim.mat')
rhos1_eq = plot([1.1*min(U_cl_vec)*U_dim 1.1*max(U_cl_vec)*U_dim], ...
                [rhos1_e_dim/rhos_0_dim rhos1_e_dim/rhos_0_dim],'--k','Linewidth',2)
rhos2_of_U = plot(U_cl_vec(1:it+1)*U_dim,rhos2_cl(1:it+1),'-','color',[0 .5 0],'LineWidth',2)
load('rhos2_e_dim.mat')
rhos2_eq = plot([1.1*min(U_cl_vec)*U_dim 1.1*max(U_cl_vec)*U_dim], ...
                [rhos2_e_dim/rhos_0_dim rhos2_e_dim/rhos_0_dim],'--','color',[0 .5 0],'Linewidth',2)
set(gca,'xlim',[1.1*min(U_cl_vec(1:it+1))*U_dim 1.1*max(U_cl_vec(1:it+1))*U_dim], ...
        'ylim',[min([rhos1_e_dim/rhos_0_dim,rhos2_e_dim/rhos_0_dim, ...
                     min(rhos1_cl(1:it+1)),min(rhos2_cl(1:it+1))])-.1 ...
                max([rhos1_e_dim/rhos_0_dim,rhos2_e_dim/rhos_0_dim, ...
                     max(rhos1_cl(1:it+1)),max(rhos2_cl(1:it+1))])+.1])
ylabel('$\ \ \ \ \frac{\rho^s_i}{\rho^s_{(0)}}$','interpreter','LaTeX','rotation',0,'FontSize',48)
grid on
legend([theta_of_U,theta_eq,rhos1_of_U,rhos1_eq,rhos2_of_U,rhos2_eq], ...
       '$\theta_c$','$\theta_{c,eq}$','$\rho^s_1$','$\rho^s_{1,eq}$', ...
       '$\rho^s_2$','$\rho^s_{1,eq}$', ...
       'interpreter','Latex','FontSize',32,'Location','East')
title('Surface densities and angle at contact line')


%Derivatives with respect to xi, used to find the tangent at the joint
phi1_xi = zeros(3,2);
xi(1) = 1;
xi(2) = -1;
for pp = 1:2
    phi1_xi(1,pp) = xi(pp)-.5;
    phi1_xi(2,pp) = -2*xi(pp);
    phi1_xi(3,pp) = xi(pp)+.5;
end

load('l_1.mat');
theta_c_mesh_vec = zeros(1,it+1);
for ii = 0:it
    load(['spine_lengths_',num2str(ii),'.mat'])
    [Nodes_rz,~] = nodes_relocator_split_v02(spine_lengths, ...
                                       spine_feet_locator,n_spines, ...
                                       n_spines_incr,n_v,alpha_s, ...
                                       n_nodes_const_per_spine, ...
                                       n_el_sing);
    r_le = Nodes_rz(l_1(1,:),1);
    z_le = Nodes_rz(l_1(1,:),2);
    tan1_r = (phi1_xi(:,2)')*r_le;
    tan1_z = (phi1_xi(:,2)')*z_le;
    theta_c_mesh_vec(ii+1) = theta_polar(-tan1_r,tan1_z);
end
figure
plot(time_vec(1:it+1)*T_dim*1E9,theta_c_vec(1:it+1)*180/pi,'linewidth',4)
hold on
grid on
plot(time_vec(1:it+1)*T_dim*1E9,theta_c_mesh_vec(1:it+1)*180/pi,'--','color',[0 1 0],'linewidth',4)
set(gca,'FontSize',16)
xlabel('$t[n$s$]$','interpreter','latex','FontSize',24)
ylabel('$\theta_c[^o]\ \ \ \ $','interpreter','LaTeX','rotation',0,'FontSize',24)
legend('$\theta_c$','$\theta_{c,t}$', ...
       'interpreter','Latex','FontSize',24,'Location','East')

figure
plot(time_vec(1:it+1)*T_dim*1E9,U_cl_vec(1:it+1)*U_dim,'linewidth',4)
hold on
grid on
set(gca,'FontSize',16)
xlabel('$t[n$s$]$','interpreter','latex','FontSize',24)
ylabel('$U_{cl}[m/s]\ \ \ \ \ \ \ \ \ \ $','interpreter','LaTeX','rotation',0,'FontSize',24)

A_vec = zeros(1,it+1);
A_vec(1) = 0;
for ii = 1:it
    load(['A_',num2str(ii),'.mat'])
    A_vec(ii+1) = A;
end
figure
plot(time_vec(1:it+1)*T_dim*1E9,A_vec,'LineWidth',2)
set(gca,'FontSize',16)
xlabel('$t[n$s$]$','interpreter','latex','FontSize',24)
title('A')
grid on
    