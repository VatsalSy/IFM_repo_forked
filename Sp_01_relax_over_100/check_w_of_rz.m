function check_w = check_w_of_rz(r,z,l,rho_l_minus2,c_theta,s_theta)

check_w =   l*r*rho_l_minus2*s_theta ...
          + l*z*rho_l_minus2*c_theta;