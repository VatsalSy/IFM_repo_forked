clear
close all
clc

load('theta_c_ini.mat')
load('hs.mat','hs')
load('T_dim.mat')
load('L_dim')
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
apex_z = zeros(1,it+1);
apex_w = zeros(1,it+1);
A = 0;
save('A_0.mat','A')
A_vec = zeros(1,it+1);
for ii = 1:it+1
    load(['rhos2_',num2str(ii-1),'.mat'])
    rhos2_cl(ii) = rhos2(1);
    load(['rhos1_',num2str(ii-1),'.mat'])
    rhos1_cl(ii) = rhos1(1);   
    load(['spine_lengths_',num2str(ii-1),'.mat'])
    apex_z(ii) = spine_lengths(end);
    load(['ws1_',num2str(ii-1),'.mat'])
    apex_w(ii) = ws1(end);
    load(['A_',num2str(ii-1),'.mat'])
    A_vec(ii) = A;
end
save('apex_z.mat','apex_z')
save('apex_z.mat','apex_z')
save('A_vec.mat','A_vec')
% load(['spine_lengths_',num2str(0),'.mat'])
% ax_symm_minus1 = spine_lengths(end);
% load(['spine_lengths_',num2str(1),'.mat'])
% ax_symm = spine_lengths(end);
% apex_w(2) = (ax_symm-ax_symm_minus1)/(time_vec(2)-time_vec(1));
% for ii = 2:it
%     ax_symm_minus2 = ax_symm_minus1;
%     ax_symm_minus1 = ax_symm;
%     load(['spine_lengths_',num2str(ii),'.mat'])
%     ax_symm = spine_lengths(end);
%     delta_t = time_vec(it+1)-time_vec(it);
%     delta_t_prev = time_vec(it)-time_vec(it-1);
%     ratio_n = delta_t/delta_t_prev;
%     a_n = (1+2*ratio_n)/(1+ratio_n);
%     a_n_minus1 = -(1+ratio_n);
%     a_n_minus2 = ratio_n^2/(1+ratio_n);
%     apex_w(ii+1) = (a_n*ax_symm+a_n_minus1*ax_symm_minus1+a_n_minus2*ax_symm_minus2)/delta_t;
% end
load('theta_c_vec.mat','theta_c_vec')
figure
yyaxis left
theta_of_U = plot(U_cl_vec(1:it+1),theta_c_vec(1:it+1)*180/pi,'-b','linewidth',4)
hold on
theta_eq = plot([1.1*min(U_cl_vec) 1.1*max(U_cl_vec)],[theta_c_eq theta_c_eq]*180/pi,'--b','Linewidth',2)
set(gca,'xlim',[1.1*min(U_cl_vec(1:it+1)) 1.1*max(U_cl_vec(1:it+1))], ...
        'ylim',180*[min(theta_c_eq,min(theta_c_vec(1:it+1)))-10/180 ...
                    max([theta_c_eq,theta_c_vec(1:it+1)])+10/180]/pi,'FontSize',32)
ylabel('$\theta_c [^o]\ \ \ \ $','interpreter','LaTeX','rotation',0,'FontSize',48)
xlabel('$\rho\nu U_{cl}/\sigma_{1,e}$','interpreter','latex','FontSize',40)
arrow_blue = annotation('arrow','Position',[.6 .8 .1 -.028],'color',[0 0 1],'LineWidth',4,'HeadLength',30,'HeadStyle','cback3')
arrow_green = annotation('arrow','Position',[.54 .315 .1 .012],'color',[0 .5 0],'LineWidth',4,'HeadLength',30,'HeadStyle','cback3')
arrow_black = annotation('arrow','Position',[.34 .215 -.1 .011],'color',[0 0 0],'LineWidth',4,'HeadLength',30,'HeadStyle','cback3')
% legend([],'$\theta_c$','$\theta_{c,eq}$', ...
%        'interpreter','Latex','FontSize',16,'Location','SouthEast')
yyaxis right
rhos1_of_U = plot(U_cl_vec(1:it+1),rhos1_cl(1:it+1),'-k','LineWidth',4)
load('rhos_0_dim.mat')
load('rhos1_e_dim.mat')
rhos1_eq = plot([1.1*min(U_cl_vec) 1.1*max(U_cl_vec)], ...
                [rhos1_e_dim/rhos_0_dim rhos1_e_dim/rhos_0_dim],'--k','Linewidth',2)
rhos2_of_U = plot(U_cl_vec(1:it+1),rhos2_cl(1:it+1),'-','color',[0 .5 0],'LineWidth',4)
load('rhos2_e_dim.mat')
rhos2_eq = plot([1.1*min(U_cl_vec) 1.1*max(U_cl_vec)], ...
                [rhos2_e_dim/rhos_0_dim rhos2_e_dim/rhos_0_dim],'--','color',[0 .5 0],'Linewidth',2)
set(gca,'xlim',[1.1*min(U_cl_vec(1:it+1)) 1.1*max(U_cl_vec(1:it+1))], ...
        'ylim',[min([rhos1_e_dim/rhos_0_dim,rhos2_e_dim/rhos_0_dim, ...
                     min(rhos1_cl(1:it+1)),min(rhos2_cl(1:it+1))])-.1 ...
                max([rhos1_e_dim/rhos_0_dim,rhos2_e_dim/rhos_0_dim, ...
                     max(rhos1_cl(1:it+1)),max(rhos2_cl(1:it+1))])+.1], ...
        'Ycolor',[0 .5 0])
ylabel('$\ \ \ \ \frac{\rho^s}{\rho^s_{(0)}}$','interpreter','LaTeX','rotation',0,'FontSize',48)
grid on
legend([theta_of_U,theta_eq,rhos1_of_U,rhos1_eq,rhos2_of_U,rhos2_eq], ...
       '$\theta_c$','$\theta_{c,eq}$','$\rho^s_1$','$\rho^s_{1,eq}$', ...
       '$\rho^s_2$','$\rho^s_{1,eq}$', ...
       'interpreter','Latex','FontSize',40,'Location','West')
set(gcf, 'Units', 'Normalized', 'OuterPosition', [0 0 1 1]);
% title('Surface densities and angle at contact line')
print('-depsc','summary_cl.eps')
pause(2)


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
index_cut = find(time_vec>200,1);
time_vec_aux = [time_vec(1:index_cut-1),200];
A_vec_aux = [A_vec(1:index_cut-1),(A_vec(index_cut-1)+A_vec(index_cut))/2];
figure
yyaxis left
plot(time_vec(1:it+1),theta_c_vec(1:it+1)*180/pi,'linewidth',9)
hold on
grid on
plot(time_vec(1:it+1),theta_c_mesh_vec(1:it+1)*180/pi,'-','color',[0 .8 0],'linewidth',3)
set(gca,'FontSize',16,'xlim',[0 200],'ylim',[60 160])
text(-48,165,' c)','FontSize',20)
xlabel('$ \sigma_{1,e}t/(\rho\nu R)$','interpreter','latex','FontSize',24)
ylabel('$\theta_c[^o]\ \ \ \ $','interpreter','LaTeX','rotation',90,'FontSize',24)
yyaxis right
plot(time_vec_aux,A_vec_aux,'LineWidth',2,'color',[0 0 0])
set(gca,'ylim',[-2000 4000],'Ycolor','k')
legend('$\theta_c$','$\theta_{c,m}$','$A$', ...
       'interpreter','Latex','FontSize',24,'Location','NorthEast')
axes('Position',[.3 .55 .3 .3])
box on
yyaxis left
plot(time_vec(1:it+1),theta_c_vec(1:it+1)*180/pi,'linewidth',9)
hold on
grid on
plot(time_vec(1:it+1),theta_c_mesh_vec(1:it+1)*180/pi,'-','color',[0 .8 0],'linewidth',3)
set(gca,'ylim',[86 98])
yyaxis right
plot(time_vec(1:index_cut-1),A_vec(1:index_cut-1),'LineWidth',2,'color',[0 0 0])
set(gca,'xlim',[0.08 .12]/(T_dim*1E6),'Ycolor','k')
print('-depsc','theta_c_f_of_time.eps')
pause(2)
   
figure
index_cut = find(time_vec>30,1);
time_vec_aux = [time_vec(1:index_cut-1),30];
A_vec_aux = [A_vec(1:index_cut-1),(A_vec(index_cut-1)+A_vec(index_cut))/2];
yyaxis left
plot(time_vec(1:it+1),theta_c_vec(1:it+1)*180/pi,'linewidth',9)
hold on
grid on
plot(time_vec(1:it+1),theta_c_mesh_vec(1:it+1)*180/pi,'-','color',[0 .8 0],'linewidth',3)
set(gca,'FontSize',16,'xlim',[0 30],'ylim',[60 160])
text(-8,165,' d)','FontSize',20)
xlabel('$\sigma_{1,e}t/(\rho\nu R)$','interpreter','latex','FontSize',24)
ylabel('$\ $','interpreter','LaTeX','rotation',90,'FontSize',24)
yyaxis right
plot(time_vec_aux,A_vec_aux,'LineWidth',2,'color',[0 0 0])
set(gca,'Ycolor',[0 0 0])
ylabel('$A$','interpreter','LaTeX','rotation',90,'FontSize',30)
legend('$\theta_c$','$\theta_{c,m}$','$A$', ...
       'interpreter','Latex','FontSize',24,'Location','NorthEast')
print('-depsc','theta_c_f_of_time_zoom.eps')
pause(2)

figure
yyaxis left
plot(time_vec(1:it+1),U_cl_vec(1:it+1),'linewidth',4)
hold on
grid on
set(gca,'FontSize',16,'xlim',[0 200],'ylim',[0 .12])
text(-43,.125,'a)','FontSize',20)
% xlabel('$t[n$s$]$','interpreter','latex','FontSize',24)
ylabel('$\rho\nu U_{cl}/\sigma_{1,e}$','interpreter','LaTeX','rotation',90,'FontSize',24)
yyaxis right
plot(time_vec(1:it+1),apex_w(1:it+1),'linewidth',2,'color',[0 .5 0])
% ylabel('$z_a\,[\mu m/s]$','interpreter','LaTeX','rotation',90,'FontSize',24)
set(gca,'Ycolor',[0 .5 0],'ylim',[-.1 .02])
legend('$U_{cl}$','$w^s_{1,a}$','interpreter','Latex','FontSize',24,'location','East')
print('-depsc','Ucl_f_of_time.eps')
pause(2)

figure
yyaxis left
plot(time_vec(1:it+1),U_cl_vec(1:it+1),'linewidth',4)
hold on
grid on
text(-8,.125,'b)','FontSize',20)
set(gca,'FontSize',16,'xlim',[0 30],'ylim',[0 .12])
% xlabel('$t[n$s$]$','interpreter','latex','FontSize',24)
ylabel('$\ $','interpreter','LaTeX','rotation',90,'FontSize',40)
yyaxis right
plot(time_vec(1:it+1),apex_w(1:it+1),'linewidth',2,'color',[0 .5 0])
ylabel('$\rho\nu w^s_{1,a}/\sigma_{1,e}$','interpreter','LaTeX','rotation',90,'FontSize',24)
set(gca,'Ycolor',[0 .5 0],'ylim',[-.1 .02])
legend('$U_{cl}$','$w^s_{1,a}$','interpreter','Latex','FontSize',24,'location','East')
print('-depsc','Ucl_f_of_time_zoom.eps')
pause(2)

% figure
% plot(time_vec(1:it+1)*T_dim*1E9,A_vec,'LineWidth',2)
% set(gca,'FontSize',16)
% xlabel('$t[n$s$]$','interpreter','latex','FontSize',24)
% title('A')
% grid on
    