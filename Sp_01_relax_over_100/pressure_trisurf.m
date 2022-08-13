close all
clear
clc


acute_theta_c = 125; %in degrees
obtuse_theta_c = 135; %in degrees

%%

%loading necessary variables
load('spine_feet_locator.mat')
load('n_spines.mat')
load('n_spines_incr.mat')
load('n_v.mat')
load('alpha_s.mat')
load('n_nodes_const_per_spine.mat')
load('n_el_sing.mat')
load('n_p.mat')
load('lp.mat')
load('l.mat')

%finding the iterations 
load('theta_c_vec.mat')
[~,index_acute] = min(abs(theta_c_vec - acute_theta_c*pi/180));
disp(['theta_c_acute = ',num2str(theta_c_vec(index_acute))])
[~,index_obtuse] = min(abs(theta_c_vec - obtuse_theta_c*pi/180));
disp(['theta_c_obtuse = ',num2str(theta_c_vec(index_obtuse))])


load(['spine_lengths_',num2str(index_acute-1),'.mat'])
[Nodes_rz,~] ...
        = nodes_relocator_split_v02(spine_lengths, ...
                                    spine_feet_locator,n_spines, ...
                                    n_spines_incr,n_v,alpha_s, ...
                                    n_nodes_const_per_spine, ...
                                    n_el_sing);
pNodes_rz = zeros(n_p,2);
for ss = 1:n_p
    [index_e,index_n] = find(lp == ss,1);
    pNodes_rz(ss,:) = Nodes_rz(l(index_e,index_n),:);
end
save('pNodes_rz.mat','pNodes_rz')                               
                                
                                
load(['p_',num2str(index_acute-1),'.mat'])

figure 
trisurf(lp(n_el_sing+1:40,:),pNodes_rz(:,1),pNodes_rz(:,2),p)
