
    
        %% test w_r derivatives
        delta_s = .0000001;
        r_plus  = r_pp + delta_s;
        z_plus  = z_pp + delta_s;
        r_minus = r_pp - delta_s;
        z_minus = z_pp - delta_s;
        %w_rr
        %forward perturbation
        theta_pol = theta_polar(r_plus,z_pp);
        phi_pp = pi - theta_c - theta_pol;
        theta_pp = lambda_eig*phi_pp;
        c_theta_pp = cos(theta_pp);
        s_theta_pp = sin(theta_pp);
        rho_l_minus2_pp = (r_plus^2+z_pp^2)^((lambda_eig-2)/2);
        rho_l_minus4_pp = (r_plus^2+z_pp^2)^((lambda_eig-4)/2);
        rho_l_minus6_pp = (r_plus^2+z_pp^2)^((lambda_eig-6)/2);
        ln_rho2_pp = log(r_plus^2+z_pp^2);
        w_r_plus(pp) = dr_check_w_of_rz(r_plus,z_pp,lambda_eig, ...
                                        rho_l_minus2_pp, ...
                                        rho_l_minus4_pp, ...
                                        c_theta_pp,s_theta_pp);
        %backward perturbation
        theta_pol = theta_polar(r_minus,z_pp);
        phi_pp = pi - theta_c - theta_pol;
        theta_pp = lambda_eig*phi_pp;
        c_theta_pp = cos(theta_pp);
        s_theta_pp = sin(theta_pp);
        rho_l_minus2_pp = (r_minus^2+z_pp^2)^((lambda_eig-2)/2);
        rho_l_minus4_pp = (r_minus^2+z_pp^2)^((lambda_eig-4)/2);
        rho_l_minus6_pp = (r_minus^2+z_pp^2)^((lambda_eig-6)/2);
        ln_rho2_pp = log(r_minus^2+z_pp^2);
        w_r_minus(pp) = dr_check_w_of_rz(r_minus,z_pp,lambda_eig, ...
                                         rho_l_minus2_pp, ...
                                         rho_l_minus4_pp, ...
                                         c_theta_pp,s_theta_pp);
        w_rr_num(pp) = (w_r_plus(pp)-w_r_minus(pp))/(2*delta_s);
        %w_rz
        %forward perturbation
        theta_pol = theta_polar(r_pp,z_plus);
        phi_pp = pi - theta_c - theta_pol;
        theta_pp = lambda_eig*phi_pp;
        c_theta_pp = cos(theta_pp);
        s_theta_pp = sin(theta_pp);
        rho_l_minus2_pp = (r_pp^2+z_plus^2)^((lambda_eig-2)/2);
        rho_l_minus4_pp = (r_pp^2+z_plus^2)^((lambda_eig-4)/2);
        rho_l_minus6_pp = (r_pp^2+z_plus^2)^((lambda_eig-6)/2);
        ln_rho2_pp = log(r_pp^2+z_plus^2);
        w_r_plus(pp) = dr_check_w_of_rz(r_pp,z_plus,lambda_eig, ...
                                        rho_l_minus2_pp, ...
                                        rho_l_minus4_pp, ...
                                        c_theta_pp,s_theta_pp);
        %backward perturbation
        theta_pol = theta_polar(r_pp,z_minus);
        phi_pp = pi - theta_c - theta_pol;
        theta_pp = lambda_eig*phi_pp;
        c_theta_pp = cos(theta_pp);
        s_theta_pp = sin(theta_pp);
        rho_l_minus2_pp = (r_pp^2+z_minus^2)^((lambda_eig-2)/2);
        rho_l_minus4_pp = (r_pp^2+z_minus^2)^((lambda_eig-4)/2);
        rho_l_minus6_pp = (r_pp^2+z_minus^2)^((lambda_eig-6)/2);
        ln_rho2_pp = log(r_pp^2+z_minus^2);
        w_r_minus(pp) = dr_check_w_of_rz(r_pp,z_minus,lambda_eig, ...
                                         rho_l_minus2_pp, ...
                                         rho_l_minus4_pp, ...
                                         c_theta_pp,s_theta_pp);
        w_rz_num(pp) = (w_r_plus(pp)-w_r_minus(pp))/(2*delta_s);