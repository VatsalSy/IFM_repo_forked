function [new_nodes,R_spines] ...
             = nodes_relocator_split_v02(spine_lengths, ...
                                         spine_feet_locator,n_spines, ...
                                         n_incr,n_nodes,alpha_s, ...
                                         n_nodes_const_per_spine, ...
                                         n_el_sing)

f = spine_lengths(1);
%Spine feet (determined using James' rule)
spine_feet = f*spine_feet_locator;

%Chi value for each spine
chi_spines = zeros(1,n_spines);
chi_spines(1) = Inf;
chi_spines(end) = 0;
for ii = 2:n_spines-1
    fun = @(chi) foot_of_chi(chi,f)-spine_feet(ii);
    chi_spines(ii) = fzero(fun,5);
end

R_spines = zeros(1,n_spines);
R_spines(2) = abs(f*(spine_feet_locator(1)-spine_feet_locator(2)));
R_spines(3) = abs(f*(spine_feet_locator(1)-spine_feet_locator(3)));
R_spines(end) = inf;
for ii = 4:n_spines-1
    R_spines(ii) = f/sinh(chi_spines(ii));
end

kk = 1;
new_nodes = zeros(n_nodes,2);
new_nodes(1,:) = [f,0];
new_nodes(2,1) = spine_feet(2);
new_nodes(2,2) = 0;
for kk = 3:n_el_sing+2
    %angle in the constant chi circle for each node 
    %(arc length over radius)
    theta = alpha_s(kk)*spine_lengths(2)/R_spines(2);
    %x of the node is the foot of the spine plus the circle's radius 
    %minus (radius*cos(theta))
    new_nodes(kk,1) = spine_feet(2) + R_spines(2)*(1-cos(theta));
    new_nodes(kk,2) = R_spines(2)*sin(theta);
end
for ii = 3:n_incr
    kk = kk+1;
    new_nodes(kk,1) = spine_feet(ii);
    new_nodes(kk,2) = 0;
    for jj = 2:2*n_el_sing+ii-2
        kk = kk+1;
        %angle in the constant chi circle for each node
        theta = alpha_s(kk)*spine_lengths(ii)/R_spines(ii); 
        new_nodes(kk,1) = spine_feet(ii) + R_spines(ii)*(1-cos(theta));
        new_nodes(kk,2) = R_spines(ii)*sin(theta);
    end
end
for ii = n_incr+1:n_spines-1
    kk = kk+1;
    new_nodes(kk,1) = spine_feet(ii);
    new_nodes(kk,2) = 0;
    for jj = 2:n_nodes_const_per_spine
        kk = kk+1;
        %angle in the constant chi circle for each node
        theta = alpha_s(kk)*spine_lengths(ii)/R_spines(ii); 
        new_nodes(kk,1) = spine_feet(ii) + R_spines(ii)*(1-cos(theta));
        new_nodes(kk,2) = R_spines(ii)*sin(theta);
    end
end
new_nodes(end-n_nodes_const_per_spine+1:end,1) = zeros(1,n_nodes_const_per_spine);
new_nodes(end-n_nodes_const_per_spine+1:end,2) ...
    = 0:spine_lengths(end)/(n_nodes_const_per_spine-1):spine_lengths(end);

