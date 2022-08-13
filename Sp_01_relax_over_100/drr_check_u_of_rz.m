function drr_check_u = drr_check_u_of_rz(r,z,l,rho_l_minus4, ...
                                         rho_l_minus6, ...
                                         c_theta,s_theta)
                                    
drr_check_u = -   l*(l-2)      *z    *rho_l_minus4*s_theta ...
              -   l*(l-2)*(l-4)*r*r*z*rho_l_minus6*s_theta ...
              -   l* l   *(l-2)*r*z*z*rho_l_minus6*c_theta ...
              -   l* l   *(l-4)*r*z*z*rho_l_minus6*c_theta ...
              +   l* l   * l   *z*z*z*rho_l_minus6*s_theta ...
              +   l*(l-2)      *r    *rho_l_minus4*c_theta ...
              -   l* l         *z    *rho_l_minus4*s_theta ...
              + 2*l*(l-2)      *r    *rho_l_minus4*c_theta ...
              +   l*(l-2)*(l-4)*r*r*r*rho_l_minus6*c_theta ...
              -   l* l   *(l-2)*r*r*z*rho_l_minus6*s_theta ...
              -   l* l         *z    *rho_l_minus4*s_theta ...
              -   l* l   *(l-4)*r*r*z*rho_l_minus6*s_theta ...
              -   l* l   * l   *r*z*z*rho_l_minus6*c_theta;
            
            
            
            