function check_u = check_u_of_rz(r,z,l,rho_l_minus2,c_theta,s_theta)

check_u = - l*z*rho_l_minus2*s_theta ...
          + l*r*rho_l_minus2*c_theta;