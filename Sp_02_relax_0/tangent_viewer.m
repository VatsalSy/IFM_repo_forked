close all
clear 
clc

load('it.mat')
load('spine_feet_locator.mat')
load('n_spines.mat')
load('n_spines_incr.mat')
load('n_v.mat')
load('phi1_lG.mat')
load('phi1_xi_lG.mat')
load('n_lGaussian_Q.mat')
load('n1_el.mat')
load('l_1.mat')
load('alpha_1.mat')

r_lGaussian = zeros(1,n_lGaussian_Q);
z_lGaussian = zeros(1,n_lGaussian_Q);
r_prime_lGaussian = zeros(1,n_lGaussian_Q);
z_prime_lGaussian = zeros(1,n_lGaussian_Q);
det_Jle = zeros(1,n_lGaussian_Q);

figure
hold on
for step = 0:it
    load(['spine_lengths_',num2str(step),'.mat'])
    [Nodes_rz,~] = nodes_relocator_split_v01(spine_lengths, ...
                                             spine_feet_locator, ...
                                             n_spines,n_spines_incr,n_v);
    for ee = 1:n1_el 
        r_le = Nodes_rz(l_1(ee,:),1);
        z_le = Nodes_rz(l_1(ee,:),2);
        for pp = 1:n_lGaussian_Q
            r_lGaussian(pp) = (phi1_lG(:,pp)')*r_le;
            z_lGaussian(pp) = (phi1_lG(:,pp)')*z_le;
            %Finding r' = \partial_\eta and z' = \partial_\eta z at 
            %each Gaussian quadrature point
            r_prime_lGaussian(pp) = (phi1_xi_lG(:,pp)')*r_le;
            z_prime_lGaussian(pp) = (phi1_xi_lG(:,pp)')*z_le;
            %Finding the determinant of the Jacobian of the 
            %isoparametric map at each Gaussian quadrature point
            det_Jle(pp) = sqrt(  r_prime_lGaussian(pp)^2 ...
                               + z_prime_lGaussian(pp)^2 ...
                              );
        end
        quiver(r_lGaussian,z_lGaussian, ...
               r_prime_lGaussian./det_Jle,z_prime_lGaussian./det_Jle)
        quiver(r_lGaussian,z_lGaussian, ...
               -alpha_1*z_prime_lGaussian./det_Jle, ...
               alpha_1*r_prime_lGaussian./det_Jle)
        axis equal
        pause
    end
    hold off
end
    