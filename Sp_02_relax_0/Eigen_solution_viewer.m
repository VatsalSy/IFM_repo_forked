close all
clear
clc


theta_c = 120  *pi/180;

%%
A = 1;
lambda_eig = pi/theta_c;

t_ind = 0;
r_ind = 0;
for tt = 0:pi/50:pi
    disp(tt/pi)
    t_ind = t_ind + 1;
    for rr = .1:.1:2
        r_ind = r_ind +1;
        r_pp = rr*cos(tt);
        z_pp = rr*sin(tt);
        x_mat(r_ind,t_ind) = r_pp;
        z_mat(r_ind,t_ind) = z_pp;
        theta_pol = theta_polar(r_pp,z_pp);
        phi_pp = pi - theta_c - theta_pol;
        theta_pp = lambda_eig*phi_pp;
        c_theta_pp = cos(theta_pp);
        s_theta_pp = sin(theta_pp);
        rho_l_minus2_pp = (r_pp^2+z_pp^2)^((lambda_eig-2)/2);
%         rho_l_minus4_pp = (r_pp^2+z_pp^2)^((lambda_eig-4)/2);
        u(r_ind,t_ind) = check_u_of_rz(r_pp,z_pp,lambda_eig,rho_l_minus2_pp, ...
                                       c_theta_pp,s_theta_pp);
        w(r_ind,t_ind) = check_w_of_rz(r_pp,z_pp,lambda_eig,rho_l_minus2_pp, ...
                                       c_theta_pp,s_theta_pp);
    end
end

figure
quiver(x_mat,z_mat,u,w)
hold on
plot([-2 0 -2*cos(theta_c)],[0 0 2*sin(theta_c)],'lineWidth',4)
        