function theta_pol = theta_polar(r,z)


if z > 0
    theta_pol = mod(atan(z/r),pi);
elseif z < 0 && r > 0
    theta_pol = atan(z/r);
elseif z < 0 && r < 0
    theta_pol = pi + atan(z/r);
elseif z == 0 && r > 0
    theta_pol = 0;
elseif z == 0 && r < 0
    theta_pol = pi;
elseif r == 0 && z > 0
    theta_pol = pi/2;
elseif r == 0 && z < 0
    theta_pol = -pi/2;
elseif r == 0 && z == 0
    theta_pol = 0;
end
% 
% if z <= 0 
%     if r >= 0
%         theta_pol = atan(z/r);
%     else 
%         theta_pol = pi+atan(z/r);
%     end
% else
%     theta_pol = mod(atan(z/r),pi);
% end
    