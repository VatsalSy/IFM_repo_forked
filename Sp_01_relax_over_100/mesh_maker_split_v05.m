close all
clear
clc

%Number of spines
load('n_spines.mat'); %this needs to be an odd number
load('n_spines_incr.mat') 
load('n_v.mat')
load('spine_feet_nodes.mat','spine_feet_nodes')
load('spine_tip_nodes.mat','spine_tip_nodes')
load('n_spines_r.mat','n_spines_r')
load('n_el_sing.mat')
load('n_v_incr.mat','n_v_incr')
load('n_nodes_const_per_spine.mat','n_nodes_const_per_spine')
%for plotting
load('Nodes_rz.mat')

%Calculating the resulting number of elements
n_el_incr = n_el_sing ...
            +((2*n_el_sing+1+2*n_el_sing+n_spines_incr-4)/2) ...
             *(n_spines_incr-3)/2; %Number of 
%elements in the region where the number of elements to the left of each 
%spine is one more than the number to the right of it.
n_el = n_el_incr ...
       + (n_spines-n_spines_incr)*(n_nodes_const_per_spine-1)/2; %the 
%total number of elements. 
save('n_el.mat','n_el')%is equal to the sum of nodes in every odd-numbered 
%spine but the last one

n_spines_near = n_spines_r;%This has to be an odd number and it's ussually 
%chosen to be less than n_spines_r + n_spines_mr
save('n_spines_near.mat','n_spines_near')
n_spines_far = n_spines-n_spines_near+1;
save('n_spines_far.mat','n_spines_far')

%Calulating the number of nodes in the near field (region that could
%require the inclusion and removal of the eigen solution as shown in S&S
%2011 CMAME)
if n_spines_near > n_spines_incr
    n_v_near = n_v_incr ...
               + (n_spines_near-n_spines_incr) ...
                 *n_nodes_const_per_spine;
    n_v_sep = n_nodes_const_per_spine;
    n_v_start_far = n_v_near - n_v_sep + 1;
    n_v_far = n_v - n_v_near + n_v_sep;
    %Calculating the number of elements in the near field region
    n_el_near = n_el_sing ...
                + (n_spines_incr-3)*n_el_sing ...
                + (1+n_spines_incr-4)*(n_spines_incr-3)/4 ...
                + (n_spines_near-n_spines_incr) ...
                  *(n_nodes_const_per_spine-1)/2;
else
    n_v_near = 2 + n_el_sing ...
               + (2*n_el_sing + 1 + 2*n_el_sing + n_spines_near - 2) ...
                 *(n_spines_near-2)/2;    
    n_v_sep = 2*n_el_sing + n_spines_near-2;
    n_v_start_far = n_v_near - n_v_sep + 1;%the first
    %of the nodes in the region
    %that lies far from the contac line
    n_v_far = n_v - n_v_near + n_v_sep;
    %Calculating the number of elements in the near field region
    n_el_near = n_el_sing ...
                + (n_spines_near-3)*n_el_sing ...
                + (1+n_spines_near-4)*(n_spines_near-3)/4;
end
save('n_v_near.mat','n_v_near')
save('n_v_start_far.mat','n_v_start_far')
save('n_v_far.mat','n_v_far')
save('n_el_near.mat','n_el_near')
save('n_v_sep.mat','n_v_sep')

%Here we create the l map, which maps the pair (element number,local 
% velocity node number) to its global node number
l = zeros(n_el,6);%matrix whose i-th row contains the global-node 
%numbers of the i-th element, the jj-th column corresponds to the node
%number of the jj-th local element node
S = zeros(n_el,3); %matrix listing the spines that go through each element
% each row of this array contains the global-spine number of the spines
%that contain one or more nodes found in this element, except for spine 1 
%which affects every single element and therefore does not need to be 
%recorded. The column position indicates the local-spine number within that
%element.
%We highlight that the first element needs to be treated separately, when
%using this matrix, as it really only has two spines that matter. We could
%have accomodated this using a cell array for S, but it is probably not
%worth the performance hit that would cause
el_s_loc_nodes = cell(n_el,3); %element's spine local nodes matrix: each 
%row of this cell array corresponds to one element and each column to a 
%local-spine number within that element. In each poistion there is a vector
%that contains the local-node number of the nodes in that element that are 
%contained in the spine in question (e.g. el_s_loc_nodes{ee,ss} = [1 4 3], 
%when the ss-th spine in element ee contain local nodes 1, 4 and 3 in it.
%We highlight that both local and global spine numbering increases as
%spines are farther form the contact line
ss = 1;
for ee = 1:n_el_sing
    l(ee,:) =  [n_el_sing + 3 + 2*(ee-1), ...
                1, ...
                n_el_sing + 5 + 2*(ee-1), ...
                n_el_sing + 4 + 2*(ee-1), ...
                2 + (ee-1), ...
                3 + (ee-1)];
    S(ee,:) = ss:ss+2;
    el_s_loc_nodes{ee,1} =  2     ;
    el_s_loc_nodes{ee,2} = [5 6  ];
    el_s_loc_nodes{ee,3} = [1 3 4];
end
next_foot_node = n_el_sing + 3; %this is just to initialise the current 
%foot node counter defined below
for ss = 3:2:n_spines_incr-2 % go through every odd number spine but the 
    %first and last
    current_foot_node = next_foot_node; %this is the node number of the 
    %corresponding spine foot
    middle_foot_node = current_foot_node+2*n_el_sing+ss-2; %node number 
    %of the foot of the next spine (we call it middle because it is in 
    %between of two spines that our counter actually goes to)
    next_foot_node = middle_foot_node+2*n_el_sing+ss-1; %node number 
    %of the foot of the one after next spine
    for jj = 1:2:2*n_el_sing+ss-4 %within each chosen spine, we go 
        %through every one of its corner nodes (i.e. taking every other 
        %node) but the last. Each one of these nodes is associated to two elements to 
        %its upper left except for the last which is only associated to one
        %element
        %
        %below we load the row of l that contains the element with one node 
        %on the current odd-numbered spine and three on the next 
        %odd-numbered spine. Numbering goes counter-clock wise as shown 
        %below (picture is drawn assuming the contact line is on the right)
        %     3
        %     |\
        %    4| \6
        %     |  \
        %     |   \
        %    1--5--2
        % here it's local node 2 the one that the counter is stopping at
        ee = ee + 1;
        l(ee,:) = [next_foot_node+jj-1, ...
                   current_foot_node+jj-1, ...
                   next_foot_node+jj+1, ...
                   next_foot_node+jj, ...
                   middle_foot_node+jj-1, ...
                   middle_foot_node+jj];
        %We now proceed to load the ee-th row of the S matrix
        S(ee,:) = ss:ss+2;
        el_s_loc_nodes{ee,1} =  2     ;
        el_s_loc_nodes{ee,2} = [5 6  ];
        el_s_loc_nodes{ee,3} = [1 3 4];
        %below we load the row of l that contains the element with three 
        %nodes on the current odd-numbered spine and two on the next 
        %odd-numbered spine. Numbering goes counter-clock wise as shown 
        %below
        %     3__6_2
        %      \   |
        %       \  |5 
        %       4\ | 
        %         \|   
        %          1
        ee = ee + 1;
        l(ee,:) = [current_foot_node+jj-1, ...
                   current_foot_node+jj+1, ...
                   next_foot_node+jj+1, ...
                   middle_foot_node+jj, ...
                   current_foot_node+jj ...
                   middle_foot_node+jj+1];
        S(ee,:) = ss:ss+2;
        el_s_loc_nodes{ee,1} = [1 2 5];
        el_s_loc_nodes{ee,2} = [4 6  ];
        el_s_loc_nodes{ee,3} =  3     ;
    end
    %here we load the last element that has a single node on the spine
    %     3
    %     |\
    %    4| \6
    %     |  \
    %     |   \
    %    1--5--2
    ee = ee + 1;
    jj = jj + 2;
    l(ee,:) = [next_foot_node+jj-1, ...
               current_foot_node+jj-1, ...
               next_foot_node+jj+1, ...
               next_foot_node+jj, ...
               middle_foot_node+jj-1, ...
               middle_foot_node+jj];
    %We now proceed to load the ee-th row of the S matrix
    S(ee,:) = ss:ss+2;
    el_s_loc_nodes{ee,1} =  2     ;
    el_s_loc_nodes{ee,2} = [5 6  ];
    el_s_loc_nodes{ee,3} = [1 3 4];
    
%     ee = ee+1;
%     l(ee,:) = [next_foot_node+ss-1; ...
%                current_foot_node+ss-1; ...
%                next_foot_node+ss+1; ...
%                next_foot_node+ss; ...
%                middle_foot_node+ss-1; ...
%                middle_foot_node+ss];
%     S(ee,:) = ss:ss+2;
%     el_s_loc_nodes{ee,1} =  2     ;
%     el_s_loc_nodes{ee,2} = [5 6  ];
%     el_s_loc_nodes{ee,3} = [1 3 4];
end
%now the part with the non-increasing number of elements
for ss = n_spines_incr:2:n_spines-2 % go through every odd number spine but 
    %the last
    current_foot_node = next_foot_node; %this is the node number of the 
    %corresponding spine foot
    middle_foot_node = current_foot_node + n_nodes_const_per_spine; %node 
    %number of the foot of the next spine (we call it middle because it is 
    %in between of two spines that our counter actually goes to)
    next_foot_node = middle_foot_node + n_nodes_const_per_spine; %node 
    %number of the foot of the one after next spine
    for jj = 1:2:n_nodes_const_per_spine-2 %within each chosen spine, we go through 
        %every one of its corner nodes (i.e. taking every other node) but 
        %the last. Each one of this nodes is associated to two elements to 
        %its upper left except for the last which is only associated to one 
        %element
        %
        %below we load the row of l that contains the element with one node 
        %on the current odd-numbered spine and three on the next 
        %odd-numbered spine. Numbering goes counter-clock wise as shown 
        %below (picture is drawn assuming the contact line is on the right
        %     3
        %     |\
        %    4| \6
        %     |  \
        %     |   \
        %    1--5--2
        % here it's local node 2 the one that the counter is stopping at
        ee = ee + 1;
        l(ee,:) = [next_foot_node+jj-1, ...
                   current_foot_node+jj-1, ...
                   next_foot_node+jj+1, ...
                   next_foot_node+jj, ...
                   middle_foot_node+jj-1, ...
                   middle_foot_node+jj];
        %We now proceed to load the ee-th row of the S matrix
        S(ee,:) = ss:ss+2;
        el_s_loc_nodes{ee,1} =  2     ;
        el_s_loc_nodes{ee,2} = [5 6  ];
        el_s_loc_nodes{ee,3} = [1 3 4];
        %below we load the row of l that contains the element with three 
        %nodes on the current odd-numbered spine and two on the next 
        %odd-numbered spine. Numbering goes counter-clock wise as shown
        %below
        %     3__6_2
        %      \   |
        %       \  |5 
        %       4\ | 
        %         \|   
        %          1
        ee = ee + 1;
        l(ee,:) = [current_foot_node+jj-1, ...
                   current_foot_node+jj+1, ...
                   next_foot_node+jj+1, ...
                   middle_foot_node+jj, ...
                   current_foot_node+jj ...
                   middle_foot_node+jj+1];
        S(ee,:) = ss:ss+2;
        el_s_loc_nodes{ee,1} = [1 2 5];
        el_s_loc_nodes{ee,2} = [4 6  ];
        el_s_loc_nodes{ee,3} =  3     ;
    end
end
save('l.mat','l')
save('S.mat','S')
el_s_loc_nodes{1,1} = [];
save('el_s_loc_nodes.mat','el_s_loc_nodes')

n_u_k_near = 0; 
save('n_u_k_near.mat','n_u_k_near')
unodes_nums_k_near = []; 
save('unodes_nums_k_near.mat','unodes_nums_k_near')

n_u_uk_near = n_v_near;
save('n_u_uk_near.mat','n_u_uk_near')
unodes_nums_uk_near = 1:n_v_near;
save('unodes_nums_uk_near.mat','unodes_nums_uk_near')

n_u_k_far = n_nodes_const_per_spine; 
save('n_u_k_far.mat','n_u_k_far')
unodes_nums_k_far  = n_v-n_nodes_const_per_spine+1:n_v; 
save('unodes_nums_k_far.mat','unodes_nums_k_far')

n_u_uk_far  = n_v_far-n_nodes_const_per_spine; 
save('n_u_uk_far.mat','n_u_uk_far')
unodes_nums_uk_far  = n_v_start_far:n_v-n_nodes_const_per_spine; 
save('unodes_nums_uk_far.mat','unodes_nums_uk_far')

n_w_k_near = 0;%n_spines_near-1; 
save('n_w_k_near.mat','n_w_k_near')

n_w_uk_near = n_v_near;%n_v_near-2*n_spines_near+1; 
save('n_w_uk_near.mat','n_w_uk_near')

n_w_k_far  = 0;%n_spines-n_spines_near+1; 
save('n_w_k_far.mat','n_w_k_far')

n_w_uk_far  = n_v_far; 
save('n_w_uk_far.mat','n_w_uk_far')

wnodes_nums_uk_near = 1:n_v_near;
wnodes_nums_k_near = [];

wnodes_nums_uk_far = n_v_start_far:n_v;
wnodes_nums_k_far = [];

save('wnodes_nums_uk_near.mat','wnodes_nums_uk_near')
save('wnodes_nums_uk_far.mat','wnodes_nums_uk_far')
save('wnodes_nums_k_near.mat','wnodes_nums_k_near')
save('wnodes_nums_k_far.mat','wnodes_nums_k_far')

%plot velocity nodes in the near field (black x marker), superpose red
%circles to the nodes where u-velocity is unknown and blue diamonds to the
%nodes where u-velocity in known
figure
plot(Nodes_rz(1:n_v_near,1),Nodes_rz(1:n_v_near,2),'xk')
hold on
u_unknown_near = plot(Nodes_rz(unodes_nums_uk_near,1),Nodes_rz(unodes_nums_uk_near,2),'or');
u_known_near = plot(Nodes_rz(unodes_nums_k_near,1),Nodes_rz(unodes_nums_k_near,2),'db');
legend(u_unknown_near,'Unknown','FontSize',16)
title('u-nodes known and unknown in near field')
axis equal

%plot velocity nodes in the far field (black x marker), superpose red
%circles to the nodes where u-velocity is unknown and blue diamonds to the
%nodes where u-velocity in known
figure
plot(Nodes_rz(n_v_start_far:end,1),Nodes_rz(n_v_start_far:end,2),'xk')
hold on
u_unknown_far = plot(Nodes_rz(unodes_nums_uk_far,1),Nodes_rz(unodes_nums_uk_far,2),'or');
u_known_far = plot(Nodes_rz(unodes_nums_k_far,1),Nodes_rz(unodes_nums_k_far,2),'db');
legend([u_unknown_far,u_known_far],'Unknown','Known','FontSize',16)
title('u-nodes known and unknown in far field')
axis equal

%plot velocity nodes in the near field (black x marker), superpose red
%circles to the nodes where w-velocity is unknown and blue diamonds to the
%nodes where w-velocity in known
figure
plot(Nodes_rz(1:n_v_near,1),Nodes_rz(1:n_v_near,2),'xk')
hold on
w_unknown_near = plot(Nodes_rz(wnodes_nums_uk_near,1),Nodes_rz(wnodes_nums_uk_near,2),'or');
w_known_near = plot(Nodes_rz(wnodes_nums_k_near,1),Nodes_rz(wnodes_nums_k_near,2),'db');
legend([w_unknown_near,w_known_near],'Unknown','Known','FontSize',16)
title('w-nodes known and unknown in near field')
axis equal

%plot velocity nodes in the far field (black x marker), superpose red
%circles to the nodes where w-velocity is unknown and blue diamonds to the
%nodes where w-velocity in known
figure
plot(Nodes_rz(n_v_start_far:end,1),Nodes_rz(n_v_start_far:end,2),'xk')
hold on
w_unknown_far = plot(Nodes_rz(wnodes_nums_uk_far,1),Nodes_rz(wnodes_nums_uk_far,2),'or');
w_known_far = plot(Nodes_rz(wnodes_nums_k_far,1),Nodes_rz(wnodes_nums_k_far,2),'db');
legend([u_unknown_far,u_known_far],'Unknown','Known','FontSize',16)
title('w-nodes known and unknown in far field')
axis equal

Mr_nodes_near = unodes_nums_uk_near;
save('Mr_nodes_near.mat','Mr_nodes_near')
Mr_nodes_far = unodes_nums_uk_far;
save('Mr_nodes_far.mat','Mr_nodes_far')
n_Mr_near = n_u_uk_near;
save('n_Mr_near.mat','n_Mr_near')
n_Mr_far = n_u_uk_far;
save('n_Mr_far.mat','n_Mr_far')

Mz_nodes_near = wnodes_nums_uk_near;
save('Mz_nodes_near.mat','Mz_nodes_near')
Mz_nodes_far = wnodes_nums_uk_far;
save('Mz_nodes_far.mat','Mz_nodes_far')
n_Mz_near = n_w_uk_near;
save('n_Mz_near.mat','n_Mz_near')
n_Mz_far = n_w_uk_far;
save('n_Mz_far.mat','n_Mz_far')

figure
plot(Nodes_rz(1:n_v_near,1),Nodes_rz(1:n_v_near,2),'xk')
hold on
Mr_needed_near = plot(Nodes_rz(Mr_nodes_near,1),Nodes_rz(Mr_nodes_near,2),'or');
legend(Mr_needed_near,'$\tilde{M}^r$ node','interpreter','latex','FontSize',16)
title('Mr-nodes near')
axis equal

figure
plot(Nodes_rz(n_v_start_far:n_v,1),Nodes_rz(n_v_start_far:n_v,2),'xk')
hold on
Mr_needed_far = plot(Nodes_rz(Mr_nodes_far,1),Nodes_rz(Mr_nodes_far,2),'or');
legend(Mr_needed_far,'$\tilde{M}^r$ node','interpreter','latex','FontSize',16)
title('Mr-nodes far')
axis equal

figure
plot(Nodes_rz(1:n_v_near,1),Nodes_rz(1:n_v_near,2),'xk')
hold on
Mz_needed_near = plot(Nodes_rz(Mz_nodes_near,1),Nodes_rz(Mz_nodes_near,2),'or');
legend(Mz_needed_near,'$\tilde{M}^z$ node','interpreter','latex','FontSize',16)
title('Mz-nodes near')
axis equal

figure
plot(Nodes_rz(n_v_start_far:n_v,1),Nodes_rz(n_v_start_far:n_v,2),'xk')
hold on
Mz_needed_far = plot(Nodes_rz(Mz_nodes_far,1),Nodes_rz(Mz_nodes_far,2),'or');
legend(Mz_needed_far,'$\tilde{M}^z$ node','interpreter','latex','FontSize',16)
title('Mz-nodes far')
axis equal


%lp function (map from local node number on triangular element to 
%pressure-node number)
lp = zeros(n_el,3);%matrix whose i-th row contains the pressure-node 
%numbers of the i-th element, the jj-th column corresponds to the node
%number of the jj-th local element node
for ee = 1:n_el_sing
    lp(ee,:) =  [2 + ee-1, ...
                1, ...
                3 + ee-1];
end
next_foot_node = 2; %this is just to initialise the current 
%foot node counter defined below
for ss = 3:2:n_spines_incr-2 % go through every odd number spine but the last
    current_foot_node = next_foot_node;%this is the pressure-node number of
    %the corresponding spine foot
    next_foot_node = current_foot_node+ n_el_sing + (ss-1)/2; %pressure-node number of 
    %the foot of the one after next spine
    for jj = 1:n_el_sing + (ss-3)/2 %within each chosen spine, we go through every one  
        % its corner nodes but the last. Each of these nodes is associated
        % to two elements to its upper left, but the last is only 
        %associated to one element
        %
        %first we load the row of lp that contains the element with a 
        %single node on the current odd-numbered spine and two on the next 
        %odd-numbered spine. Numbering goes counter-clock wise as shown 
        %below (picture is drawn assuming the contact line is on the right)
        %     3
        %     |\
        %     | \
        %     |  \
        %     |   \
        %    1-----2
        % here it's local node 2 the one that the jj counter is stopping at
        ee = ee + 1;
        lp(ee,:) = [next_foot_node+jj-1, ...
                    current_foot_node+jj-1, ...
                    next_foot_node+jj];
        %below we load the row of lp that contains the element with two 
        %nodes on the current odd-numbered spine and two on the next 
        %odd-numbered spine. Numbering goes counter-clock wise as shown 
        %below
        %     3____2
        %      \   |
        %       \  |
        %        \ | 
        %         \|   
        %          1
        ee = ee + 1;
        lp(ee,:) = [current_foot_node+jj-1, ...
                      current_foot_node+jj, ...
                      next_foot_node+jj];
    end
    %here we load the last element that has a single node on the spine
    %     3
    %     |\
    %     | \
    %     |  \
    %     |   \
    %    1-----2
    ee = ee + 1;
    jj = jj + 1;
    lp(ee,:) = [next_foot_node+jj-1, ...
                current_foot_node+jj-1, ...
                next_foot_node+jj];
%     ee = ee + 1;
%     lp(ee,:) = [next_foot_node+(ss-1)/2, ...
%                 current_foot_node+(ss-1)/2, ...
%                 next_foot_node+(ss-1)/2+1];
end
%doing the part where the number of elements per spine no longer increases
for ss = n_spines_incr:2:n_spines-2 % go n_nodes_incrthrough every odd 
    %number spine but the last
    current_foot_node = next_foot_node;%this is the pressure-node number of
    %the corresponding spine foot
    next_foot_node = current_foot_node ...
                     + (n_nodes_const_per_spine+1)/2; %pressure-node 
    %number of the foot of the one after next spine
    for jj = 1:(n_nodes_const_per_spine-1)/2 %within each chosen spine, we go through
        %every one its corner nodes but the last. Each of these nodes is 
        %associated to two elements to its upper left, but the last is only 
        %associated to one element
        %
        %first we load the row of lp that contains the element with a 
        %single node on the current odd-numbered spine and two on the next 
        %odd-numbered spine. Numbering goes counter-clock wise as shown 
        %below (picture is drawn assuming the contact line is on the right)
        %     3
        %     |\
        %     | \
        %     |  \
        %     |   \
        %    1-----2
        % here it's local node 2 the one that the jj counter is stopping at
        ee = ee + 1;
        lp(ee,:) = [next_foot_node+jj-1, ...
                    current_foot_node+jj-1, ...
                    next_foot_node+jj];
        %below we load the row of lp that contains the element with two 
        %nodes on the current odd-numbered spine and two on the next 
        %odd-numbered spine. Numbering goes counter-clock wise as shown 
        %below
        %     3____2
        %      \   |
        %       \  |
        %        \ | 
        %         \|   
        %          1
        ee = ee + 1;
        lp(ee,:) = [current_foot_node+jj-1, ...
                    current_foot_node+jj, ...
                    next_foot_node+jj];
    end
end
save('lp.mat','lp')

%number of pressure nodes
n_p = n_el_sing + 2 ...
      + (n_el_sing+1)*(n_spines_incr-3)/2 ...
      + (1+(n_spines_incr-3)/2)*(n_spines_incr-3)/4 ...
      + (n_nodes_const_per_spine+1) ...
        *(n_spines-n_spines_incr)/4;
save('n_p.mat','n_p')
if n_spines_near > n_spines_incr
    n_p_near = n_el_sing + 2 ...
               + (n_el_sing+1)*(n_spines_incr-3)/2 ...
               + (1+(n_spines_incr-3)/2)*(n_spines_incr-3)/4 ...
               + (n_nodes_const_per_spine+1) ...
                 *(n_spines_near-n_spines_incr)/4;
    n_p_sep = (n_nodes_const_per_spine+1)/2;
    n_p_start_far = n_p_near - n_p_sep + 1;
    n_p_far = n_p - n_p_near + n_p_sep;
    n_p_uk_far = n_p-n_p_near + n_p_sep;
    pnodes_nums_uk_far = n_p_start_far:n_p;
else
    n_p_near = n_el_sing + 2 ...
               + (n_el_sing+1)*(n_spines_near-3)/2 ...
               + (1+(n_spines_near-3)/2)*(n_spines_near-3)/4;
    n_p_sep = n_el_sing + (n_spines_near-1)/2;
    n_p_start_far  = n_p_near - n_p_sep + 1;
    n_p_far = n_p - n_p_near + n_p_sep;
    n_p_uk_far = n_p-n_p_near + n_p_sep;
    pnodes_nums_uk_far = n_p_start_far:n_p;
end
save('n_p_near.mat','n_p_near')
save('n_p_sep.mat','n_p_sep')
save('n_p_start_far.mat','n_p_start_far')
save('n_p_far.mat','n_p_far')
save('n_p_uk_far.mat','n_p_uk_far')
save('pnodes_nums_uk_far.mat','pnodes_nums_uk_far')
%number of pressure nodes where the pressure is unknown (i.e. everywhere
%but the free surface)

n_p_k_near = 0;
save('n_p_k_near.mat','n_p_k_near')
pnodes_nums_k_near = [];
save('pnodes_nums_k_near.mat','pnodes_nums_k_near')

n_p_uk_near = n_p_near;
save('n_p_uk_near.mat','n_p_uk_near')
pnodes_nums_uk_near = 1:n_p_near;
save('pnodes_nums_uk_near.mat','pnodes_nums_uk_near')

n_p_k_far = 0;
save('n_p_k_far.mat','n_p_k_far')

pnodes_nums_k_far = [];
save('pnodes_nums_k_far.mat','pnodes_nums_k_far')

n_C_near = n_p_uk_near;
save('n_C_near.mat','n_C_near')
C_nodes_near = pnodes_nums_uk_near;
save('C_nodes_near.mat','C_nodes_near')
n_C_far = n_p_uk_far;
save('n_C_far.mat','n_C_far')
C_nodes_far = pnodes_nums_uk_far;
save('C_nodes_far.mat','C_nodes_far')

pNodes_rz = zeros(n_p,2);
for ss = 1:n_p
    [index_e,index_n] = find(lp == ss,1);
    pNodes_rz(ss,:) = Nodes_rz(l(index_e,index_n),:);
end
save('pNodes_rz.mat','pNodes_rz')

figure
plot(Nodes_rz(1:n_v_near,1),Nodes_rz(1:n_v_near,2),'xk')
hold on
p_unknown_near = plot(pNodes_rz(pnodes_nums_uk_near,1),pNodes_rz(pnodes_nums_uk_near,2),'or');
p_known_near = plot(pNodes_rz(pnodes_nums_k_near,1),pNodes_rz(pnodes_nums_k_near,2),'db');
legend(p_unknown_near,'$C$ node','interpreter','latex','FontSize',16)
title('p-nodes known and unknown')
axis equal

figure
plot(Nodes_rz(n_v_start_far:n_v,1),Nodes_rz(n_v_start_far:n_v,2),'xk')
hold on
p_unknown_far = plot(pNodes_rz(pnodes_nums_uk_far,1),pNodes_rz(pnodes_nums_uk_far,2),'or');
p_known_far = plot(pNodes_rz(pnodes_nums_k_far,1),pNodes_rz(pnodes_nums_k_far,2),'db');
legend(p_unknown_far,'$C$ node','interpreter','latex','FontSize',16)
title('p-nodes known and unknown')
axis equal

%Boundary_1 (free surface) 
n1_el = (n_spines-1)/2;%there is one line element for every two spines but
%the last
save('n1_el.mat','n1_el')%number of line elements on boundary 1
n1_el_near = (n_spines_near-1)/2; %Number of line elements on boundary 1
%within the near zone
save('n1_el_near.mat','n1_el_near')
l_1 = zeros(n1_el,3);%the l_1 matrix stores in its k-th row the global node 
%numbers of the nodes in the k-th line element of boundary 1, column 
%numbers correspond to local line-element node number. Notice that the
%numbering convention moves away from the contact line and towards the apex
%both in the local numbering and in the boundary-1 numbering
kk = 1;
l_1(kk,:) = [1, n_el_sing+2, 3*n_el_sing+3];
last = 3*n_el_sing+3;
for ss = 3:2:n_spines_incr-2
    kk = kk+1;
    l_1(kk,1) = last;
    l_1(kk,2) = last+2*n_el_sing+ss-1;
    l_1(kk,3) = l_1(kk,2)+2*n_el_sing+ss;
    last = l_1(kk,3);
end
for ss = n_spines_incr:2:n_spines-2
    kk = kk+1;
    l_1(kk,1) = last;
    l_1(kk,2) = last+n_nodes_const_per_spine;
    l_1(kk,3) = l_1(kk,2)+n_nodes_const_per_spine;
    last = l_1(kk,3);
end
save('l_1.mat','l_1')
l1_1 = zeros(n1_el,3);%the l1_1 matrix stores on its k-th row the boundary 
%node-number the k-th line element of boundary 1, column numbers correspond
%to local line-element node number
S_1 = zeros(n1_el,3);
el1_s_loc_nodes = zeros(n1_el,3);
kk = 0;
last = 1;
for ss = 1:2:n_spines-2
    kk = kk+1;
    l1_1(kk,1) = last;
    l1_1(kk,2) = last+1;
    l1_1(kk,3) = last+2;
    S_1(kk,:) = last:last+2;
    el1_s_loc_nodes(kk,1) = 1;
    el1_s_loc_nodes(kk,2) = 2;
    el1_s_loc_nodes(kk,3) = 3;
    last = last+2;
end
save('l1_1.mat','l1_1')
save('S_1.mat','S_1')
save('el1_s_loc_nodes.mat','el1_s_loc_nodes')

%For this particular problem the exact same method goes for S_2
S_2 = S_1;
save('S_2.mat','S_2')
el2_s_loc_nodes = el1_s_loc_nodes;
save('el2_s_loc_nodes.mat','el2_s_loc_nodes')

%Boundary_2 (solid surface) stores the nodes in each element of the solid
%surface. Notice that here also the local and global numbering both
%increase as the nodes are farther away from the contact line
n2_el = n1_el;
save('n2_el.mat','n2_el')
n2_el_near = (n_spines_near-1)/2; %Number of line elements on boundary 2
%within the near zone
save('n2_el_near.mat','n2_el_near')
l_2 = zeros(n2_el,3);
kk = 1;
l_2(kk,:) = [1, 2, n_el_sing+3];
last = n_el_sing+3;
for ss = 3:2:n_spines-2
    kk = kk+1;
    %numbering from right to left
    l_2(kk,1) = last;
    l_2(kk,2) = l_1(kk,1)+1;
    l_2(kk,3) = l_1(kk,2)+1;
    last = l_2(kk,3);
end
save('l_2.mat','l_2')
l2_2 = zeros(n2_el,3);
kk = 0;
last = 1;
for ss = 1:2:n_spines-2
    kk = kk+1;
    %numbering from right to left
    l2_2(kk,1) = last;
    l2_2(kk,2) = last+1;
    l2_2(kk,3) = last+2;
    last = last+2;
end
save('l2_2.mat','l2_2')

%Boundary_3 (fluid flow) stores the nodes in each line element of the fluid
%flow boundary. Notice that here both the local numbering and the
%boundary-3 numbering increase as the nodes are farther from the solid
%interface and closer to the apex
n3_el = (n_nodes_const_per_spine-1)/2;
save('n3_el.mat','n3_el')
l_3 = zeros(n3_el,3);
kk = 0;
last = spine_feet_nodes(end);
for ss = 1:2:n_nodes_const_per_spine-2
    kk = kk+1;
    %numbering is done from bottom to top
    l_3(kk,1) = last;
    last = last+1;
    l_3(kk,2) = last;
    last = last+1;
    l_3(kk,3) = last;
end
save('l_3.mat','l_3')
l3_3 = zeros(n3_el,3);
kk = 0;
last = 1;
for ss = 1:2:n_spines-2
    kk = kk+1;
    %numbering is done from bottom to top
    l3_3(kk,1) = last;
    last = last + 1;
    l3_3(kk,2) = last;
    last = last + 1;
    l3_3(kk,3) = last;
end
save('l3_3.mat','l3_3')
S_3 = n_spines*ones((n_nodes_const_per_spine-1)/2,1);
save('S_3.mat','S_3')
el3_s_loc_nodes = cell((n_nodes_const_per_spine-1)/2,1);
for ee = 1:(n_nodes_const_per_spine-1)/2
    el3_s_loc_nodes{ee,1} = [1 2 3];
end
save('el3_s_loc_nodes.mat','el3_s_loc_nodes')


%Boundary_4 (separatix) stores the nodes in each line element of the domain
%separatrix. Notice that here both the local numbering and the
%boundary-4 numbering increase as the nodes are farther from the solid
%interface and closer to the free surface
if n_spines_near > n_spines_incr
    n4_el = (n_nodes_const_per_spine-1)/2;
    l_4 = zeros(n4_el,3);
    kk = 0;
    last = spine_feet_nodes(n_spines_near);
    for ss = 1:2:(n_nodes_const_per_spine-2)
        kk = kk+1;
        %numbering is done from bottom to top
        l_4(kk,1) = last;
        last = last + 1;
        l_4(kk,2) = last;
        last = last + 1;
        l_4(kk,3) = last;
    end
    l4_4 = zeros(n4_el,3);
    S_4 = zeros(n4_el,1);
    kk = 0;
    last = 1;
    el4_s_loc_nodes = cell(n4_el,1);
    for ss = 1:2:(n_nodes_const_per_spine-2)
        kk = kk+1;
        %numbering is done from bottom to top
        l4_4(kk,1) = last;
        last = last + 1;
        l4_4(kk,2) = last;
        last = last + 1;
        l4_4(kk,3) = last;
        S_4(kk,:) = n_spines_near;
        el4_s_loc_nodes{kk,1} = [1 2 3];
    end
else
    n4_el = (2*n_el_sing+n_spines_near-3)/2;
    l_4 = zeros(n4_el,3);
    kk = 0;
    last = spine_feet_nodes(n_spines_near);
    for ss = 1:2:(2*n_el_sing+n_spines_near-4)
        kk = kk+1;
        %numbering is done from bottom to top
        l_4(kk,1) = last;
        last = last + 1;
        l_4(kk,2) = last;
        last = last + 1;
        l_4(kk,3) = last;
    end
    l4_4 = zeros(n4_el,3);
    S_4 = zeros(n4_el,1);
    kk = 0;
    last = 1;
    el4_s_loc_nodes = cell(n4_el,1);
    for ss = 1:2:(2*n_el_sing+n_spines_near-4)
        kk = kk+1;
        %numbering is done from bottom to top
        l4_4(kk,1) = last;
        last = last + 1;
        l4_4(kk,2) = last;
        last = last + 1;
        l4_4(kk,3) = last;
        S_4(kk,:) = n_spines_near;
        el4_s_loc_nodes{kk,1} = [1 2 3];
    end
end
save('n4_el.mat','n4_el')
save('l_4.mat','l_4')
save('l4_4.mat','l4_4')
save('S_4.mat','S_4')
save('el4_s_loc_nodes.mat','el4_s_loc_nodes')

%Boundary_5 (separatix) stores the nodes in each line element of the domain
%separatrix. Notice that here both the local numbering and the
%boundary_5 numbering increase as the nodes are farther from the solid
%interface and closer to the free surface
if n_spines_near > n_spines_incr
    n5_el = n4_el;
    l_5 = zeros(n5_el,3);
    kk = 0;
    last = n_v_near;
    for ss = 1:n5_el
        kk = kk+1;
        %numbering is done from bottom to top
        l_5(kk,1) = last;
        last = last - 1;
        l_5(kk,2) = last;
        last = last - 1;
        l_5(kk,3) = last;
    end
    l5_5 = zeros(n5_el,3);
    S_5 = zeros(n5_el,1);
    el5_s_loc_nodes = cell(n4_el,1);
    kk = 0;
    last = n_v_sep;
    for ss = 1:n5_el
        kk = kk+1;
        %numbering is done from bottom to top
        l5_5(kk,1) = last;
        last = last - 1;
        l5_5(kk,2) = last;
        last = last - 1;
        l5_5(kk,3) = last;
        S_5(kk,:) = n_spines_near;
        el5_s_loc_nodes{kk,1} = [1 2 3];
    end
else
    n5_el = n4_el;
    l_5 = zeros(n5_el,3);
    kk = 0;
    last = n_v_near;
    for ss = 1:n5_el
        kk = kk+1;
        %numbering is done from bottom to top
        l_5(kk,1) = last;
        last = last - 1;
        l_5(kk,2) = last;
        last = last - 1;
        l_5(kk,3) = last;
    end
    l5_5 = zeros(n5_el,3);
    S_5 = zeros(n5_el,1);
    el5_s_loc_nodes = cell(n4_el,1);
    kk = 0;
    last = n_v_sep;
    for ss = 1:n5_el
        kk = kk+1;
        %numbering is done from bottom to top
        l5_5(kk,1) = last;
        last = last - 1;
        l5_5(kk,2) = last;
        last = last - 1;
        l5_5(kk,3) = last;
        S_5(kk,:) = n_spines_near;
        el5_s_loc_nodes{kk,1} = [1 2 3];
    end
end
save('n5_el.mat','n5_el')
save('l_5.mat','l_5')
save('l5_5.mat','l5_5')
save('S_5.mat','S_5')
save('el5_s_loc_nodes.mat','el5_s_loc_nodes')

contact_line = 1;
save('contact_line.mat','contact_line')

apex = n_v;
save('apex.mat','apex')

%Matrix T (required to calculate the Jacobian of the global to local 
%coordinates)
%Gaussian quadrature points
V1 = [-1  1];
V2 = [-1 -1];
V3 = [ 1 -1];

n_Gaussian_Q = 16; save('n_Gaussian_Q.mat','n_Gaussian_Q');
xi_eta_G = zeros(16,2);
W_G = zeros(1,16);

xi_eta_G(1,:) = V1/3+V2/3+V3/3;
W_G(1) = .1443156076777871682510911104890646;

a = [.1705693077517602066222935014914645, ...
     .0505472283170309754584235505965989, ...
     .4592925882927231560288155144941693];
p = [.1032173705347182502817915502921290, ...
     .0324584976231980803109259283417806, ...
     .0950916342672846247938961043885843];
kk = 1;
for ii = 1:3
    W_G(kk+1:kk+3) = p(ii)*ones(1,3);
    kk = kk+1;
    xi_eta_G(kk,:) =       a(ii)*V1 +       a(ii)*V2 + (1-2*a(ii))*V3;
    kk = kk+1;
    xi_eta_G(kk,:) =       a(ii)*V1 + (1-2*a(ii))*V2 +       a(ii)*V3;
    kk = kk+1;
    xi_eta_G(kk,:) = (1-2*a(ii))*V1 +       a(ii)*V2 +       a(ii)*V3;
end

b = .2631128296346381134217857862846436;
c = .0083947774099576053372138345392944;

W_G(kk+1:kk+6) = .0272303141744349942648446900739089*ones(1,6);
kk = kk+1;
xi_eta_G(kk,:) =       b*V1 +       c*V2 + (1-b-c)*V3;
kk = kk+1;
xi_eta_G(kk,:) =       b*V1 + (1-b-c)*V2 +       c*V3;
kk = kk+1;
xi_eta_G(kk,:) =       c*V1 +       b*V2 + (1-b-c)*V3;
kk = kk+1;
xi_eta_G(kk,:) =       c*V1 + (1-b-c)*V2 +       b*V3;
kk = kk+1;
xi_eta_G(kk,:) = (1-b-c)*V1 +       b*V2 +       c*V3;
kk = kk+1;
xi_eta_G(kk,:) = (1-b-c)*V1 +       c*V2 +       b*V3;

xi_G = xi_eta_G(:,1)';
save('xi_G.mat','xi_G')
eta_G = xi_eta_G(:,2)';
save('eta_G.mat','eta_G')
W_G = 2*W_G; %the weights are for a unit area triangle, so they need to be 
%multiplied by the are of the triangle
save('W_G.mat','W_G')


% n_Gaussian_Q = 9;
% save('n_Gaussian_Q.mat','n_Gaussian_Q');
% xi_G = zeros(1,n_Gaussian_Q);
% xi_G(1) = 0; xi_G(2) = 0; xi_G(3) = 0; ...
% xi_G(4) = 0.774596669241483; xi_G(5) = xi_G(4); xi_G(6) = xi_G(4); ...
% xi_G(7) = -xi_G(4); xi_G(8) = -xi_G(4); xi_G(9) = -xi_G(4);
% eta_G = zeros(1,n_Gaussian_Q);
% eta_G(1) = -0.887298334620741; eta_G(2) = -.5; eta_G(3) = -0.112701665379258;
% eta_G(4) = -0.974596669241483; eta_G(5) = eta_G(1); eta_G(6) = -0.8; 
% eta_G(7) = eta_G(6); eta_G(8) = eta_G(3); eta_G(9) = 0.574596669241483;
% W_G = zeros(1,n_Gaussian_Q);
% W_G(1) = .246913580246913; W_G(2) = .395061728395061; W_G(3) = W_G(1); 
% W_G(4) = .034784464623227; W_G(5) = .055655143397164; W_G(6) = W_G(4);
% W_G(7) = .273857510685414; W_G(8) = .438172017096662; W_G(9) = W_G(7);
% save('W_G.mat','W_G')
T = cell(1,n_Gaussian_Q);
phi_xi = zeros(6,n_Gaussian_Q);
phi_eta = zeros(6,n_Gaussian_Q);
for pp = 1:n_Gaussian_Q
    phi_xi(1,pp) = 0; 
    phi_xi(2,pp) = .5+xi_G(pp)+eta_G(pp); 
    phi_xi(3,pp) = .5+xi_G(pp);
    phi_xi(4,pp) = 1+eta_G(pp); 
    phi_xi(5,pp) = -1-eta_G(pp); 
    phi_xi(6,pp) = -1-eta_G(pp)-2*xi_G(pp);
    phi_eta(1,pp) = .5+eta_G(pp); 
    phi_eta(2,pp) = .5+xi_G(pp)+eta_G(pp); 
    phi_eta(3,pp) = 0;
    phi_eta(4,pp) = 1+xi_G(pp); 
    phi_eta(5,pp) = -1-xi_G(pp)-2*eta_G(pp); 
    phi_eta(6,pp) = -1-xi_G(pp);
    for ii = 1:6
        for jj = 1:6
            T{pp}(ii,jj) = phi_xi(ii,pp)*phi_eta(jj,pp) ...
                           -phi_eta(ii,pp)*phi_xi(jj,pp);
        end
    end
end
save('T.mat','T')

%Matrix psi_G stores the values of psi_i in each of the Guassian 
%integration points
psi_G = zeros(3,n_Gaussian_Q);
for pp = 1:n_Gaussian_Q
    psi_G(1,pp) = (1+eta_G(pp))/2;
end
for pp = 1:n_Gaussian_Q
    psi_G(2,pp) = -(xi_G(pp)+eta_G(pp))/2;
end
for pp = 1:n_Gaussian_Q
    psi_G(3,pp) = (1+xi_G(pp))/2;
end
save('psi_G.mat','psi_G')

%Matrix phi_G stores the values of phi_i in each of the Guassian 
%integration points
phi_G = zeros(6,n_Gaussian_Q);
for jj = 1:3
    for pp = 1:n_Gaussian_Q
        phi_G(jj,pp) = psi_G(jj,pp)*(2*psi_G(jj,pp)-1);
    end
end
for pp = 1:n_Gaussian_Q
    phi_G(4,pp) = 4*psi_G(1,pp)*psi_G(3,pp);
end
for jj = 5:6
    for pp = 1:n_Gaussian_Q
        phi_G(jj,pp) = 4*psi_G(jj-3,pp)*psi_G(jj-4,pp);
    end
end
save('phi_G.mat','phi_G')

%Matrices needed for line-element integration
% n_lGaussian_Q = 4;
xi_lG = [-0.9782286581460569928039, -0.8870625997680952990752, ...
         -0.7301520055740493240934, -0.5190961292068118159257, ...
         -0.2695431559523449723315, 0, 0.269543155952344972332, ...
         0.5190961292068118159257, 0.7301520055740493240934 ...
         0.887062599768095299075, 0.9782286581460569928039];
eta_lG = xi_lG;
n_lGaussian_Q = length(xi_lG);
save('n_lGaussian_Q.mat','n_lGaussian_Q')
% xi_lG = zeros(1,4);
% xi_lG(1) = -0.861136311594052; xi_lG(2) = -xi_lG(1);
% xi_lG(3) = -0.339981043584856; xi_lG(4) = -xi_lG(3);
save('xi_lG.mat','xi_lG')
% W_lG = zeros(1,4);
% W_lG(1) = 0.347854845137453; W_lG(2) = W_lG(1);
% W_lG(3) = 0.652145154862546; W_lG(4) = W_lG(3);
W_lG = [0.0556685671161736664828, 0.1255803694649046246347, ...
        0.1862902109277342514261, 0.2331937645919904799185, ...
        0.2628045445102466621807, 0.2729250867779006307145, ...
        0.262804544510246662181, 0.2331937645919904799185, ...
        0.1862902109277342514261, 0.1255803694649046246347, ...
        0.055668567116173666483];
save('W_lG.mat','W_lG')
%assesing the line-element hat functions at the Gaussian quadrature points
%for boundary 1
phi1_lG = zeros(3,n_lGaussian_Q);
for pp = 1:n_lGaussian_Q
    phi1_lG(1,pp) = xi_lG(pp)^2/2-xi_lG(pp)/2;
    phi1_lG(2,pp) = 1-xi_lG(pp)^2;
    phi1_lG(3,pp) = xi_lG(pp)^2/2+xi_lG(pp)/2;
end
save('phi1_lG.mat','phi1_lG')
%assessing the derivatives of the hat functions with respect to the one
%dimensional variable at the Gaussian quadrature points in the line
%elements of boundary1
phi1_xi_lG = zeros(3,n_lGaussian_Q);
for pp = 1:n_lGaussian_Q
    phi1_xi_lG(1,pp) = xi_lG(pp)-.5;
    phi1_xi_lG(2,pp) = -2*xi_lG(pp);
    phi1_xi_lG(3,pp) = xi_lG(pp)+.5;
end
save('phi1_xi_lG.mat','phi1_xi_lG')
phi1_xixi_lG = zeros(3,n_lGaussian_Q);
for pp = 1:n_lGaussian_Q
    phi1_xixi_lG(1,pp) = 1;
    phi1_xixi_lG(2,pp) = -2;
    phi1_xixi_lG(3,pp) = 1;
end
save('phi1_xixi_lG.mat','phi1_xixi_lG')

%assesing the line-element hat functions at the Gaussian quadrature points
%for boundary 2
phi2_lG = zeros(3,n_lGaussian_Q);
for pp = 1:n_lGaussian_Q
    phi2_lG(1,pp) = xi_lG(pp)^2/2-xi_lG(pp)/2;
    phi2_lG(2,pp) = 1-xi_lG(pp)^2;
    phi2_lG(3,pp) = xi_lG(pp)^2/2+xi_lG(pp)/2;
end
save('phi2_lG.mat','phi2_lG')
%assessing the derivatives of the hat functions with respect to the one
%dimensional variable at the Gaussian quadrature points in the line
%elements of boundary 2
phi2_eta_lG = zeros(3,n_lGaussian_Q);
for pp = 1:n_lGaussian_Q
    phi2_eta_lG(1,pp) = xi_lG(pp)-.5;%using xi here because the Gaussian 
    %int points are the same
    phi2_eta_lG(2,pp) = -2*xi_lG(pp);
    phi2_eta_lG(3,pp) = xi_lG(pp)+.5;
end
save('phi2_eta_lG.mat','phi2_eta_lG')

%assesing the line-element hat functions at the Gaussian quadrature points
%for boundary 3
phi3_lG = zeros(3,n_lGaussian_Q);
for pp = 1:n_lGaussian_Q
    phi3_lG(1,pp) = xi_lG(pp)^2/2-xi_lG(pp)/2;
    phi3_lG(2,pp) = 1-xi_lG(pp)^2;
    phi3_lG(3,pp) = xi_lG(pp)^2/2+xi_lG(pp)/2;
end
save('phi3_lG.mat','phi3_lG')
%assessing the derivatives of the hat functions with respect to the one
%dimensional variable at the Gaussian quadrature points in the line
%elements of boundary 3
phi3_xi_lG = zeros(3,n_lGaussian_Q);
for pp = 1:n_lGaussian_Q
    phi3_xi_lG(1,pp) = xi_lG(pp)-.5;%using xi here because the Gaussian 
    %int points are the same
    phi3_xi_lG(2,pp) = -2*xi_lG(pp);
    phi3_xi_lG(3,pp) = xi_lG(pp)+.5;                            
end
save('phi3_xi_lG.mat','phi3_xi_lG')

%assesing the line-element hat functions at the Gaussian quadrature points
%for boundary 4
phi4_lG = zeros(3,n_lGaussian_Q);
for pp = 1:n_lGaussian_Q
    phi4_lG(1,pp) = xi_lG(pp)^2/2-xi_lG(pp)/2;
    phi4_lG(2,pp) = 1-xi_lG(pp)^2;
    phi4_lG(3,pp) = xi_lG(pp)^2/2+xi_lG(pp)/2;
end
save('phi4_lG.mat','phi4_lG')
%assessing the derivatives of the hat functions with respect to the one
%dimensional variable at the Gaussian quadrature points in the line
%elements of boundary 3
phi4_xi_lG = zeros(3,n_lGaussian_Q);
for pp = 1:n_lGaussian_Q
    phi4_xi_lG(1,pp) = xi_lG(pp)-.5;%using xi here because the Gaussian 
    %int points are the same
    phi4_xi_lG(2,pp) = -2*xi_lG(pp);
    phi4_xi_lG(3,pp) = xi_lG(pp)+.5;                            
end
save('phi4_xi_lG.mat','phi4_xi_lG')

%assesing the line-element hat functions at the Gaussian quadrature points
%for boundary 4
phi5_lG = zeros(3,n_lGaussian_Q);
for pp = 1:n_lGaussian_Q
    phi5_lG(1,pp) = eta_lG(pp)^2/2+eta_lG(pp)/2;
    phi5_lG(2,pp) = 1-eta_lG(pp)^2;
    phi5_lG(3,pp) = eta_lG(pp)^2/2-eta_lG(pp)/2;
end
save('phi5_lG.mat','phi5_lG')
%assessing the derivatives of the hat functions with respect to the one
%dimensional variable at the Gaussian quadrature points in the line
%elements of boundary 3
phi5_eta_lG = zeros(3,n_lGaussian_Q);
for pp = 1:n_lGaussian_Q
    phi5_eta_lG(1,pp) = eta_lG(pp)+.5;%using xi here because the Gaussian 
    %int points are the same
    phi5_eta_lG(2,pp) = -2*eta_lG(pp);
    phi5_eta_lG(3,pp) = eta_lG(pp)-.5;                            
end
save('phi5_eta_lG.mat','phi5_eta_lG')


alpha_1 = 1;
save('alpha_1.mat','alpha_1')
alpha_2 = -1;
save('alpha_2.mat','alpha_2')
alpha_3 = -1;
save('alpha_3.mat','alpha_3')
alpha_4 = -1;
save('alpha_4.mat','alpha_4')
alpha_5 = -1;
save('alpha_5.mat','alpha_5')

%Creating matrix Tmat1
Tmat1 = cell(1,n_lGaussian_Q);
phi_xi = zeros(6,n_lGaussian_Q);
phi_eta = zeros(6,n_lGaussian_Q);
for pp = 1:n_lGaussian_Q
    phi_xi(1,pp) = 0; 
    phi_xi(2,pp) = -.5+xi_lG(pp); 
    phi_xi(3,pp) = .5+xi_lG(pp);
    phi_xi(4,pp) = 0; 
    phi_xi(5,pp) = 0; 
    phi_xi(6,pp) = -2*xi_lG(pp);
    phi_eta(1,pp) = -.5; 
    phi_eta(2,pp) = -.5+xi_lG(pp); 
    phi_eta(3,pp) = 0;
    phi_eta(4,pp) = 1+xi_lG(pp); 
    phi_eta(5,pp) = 1-xi_lG(pp); 
    phi_eta(6,pp) = -1-xi_lG(pp);
    for ii = 1:6
        for jj = 1:6
            Tmat1{pp}(ii,jj) = phi_xi(ii,pp)*phi_eta(jj,pp) ...
                              -phi_eta(ii,pp)*phi_xi(jj,pp);
        end
    end
end
save('Tmat1.mat','Tmat1')

%Creating matrix Tmat2
Tmat2 = cell(1,n_lGaussian_Q);
phi_xi = zeros(6,n_lGaussian_Q);
phi_eta = zeros(6,n_lGaussian_Q);
for pp = 1:n_lGaussian_Q
    phi_xi(1,pp) = 0; 
    phi_xi(2,pp) = -.5+eta_lG(pp); 
    phi_xi(3,pp) = -.5;
    phi_xi(4,pp) = 1+eta_lG(pp); 
    phi_xi(5,pp) = -1-eta_lG(pp); 
    phi_xi(6,pp) = 1-eta_lG(pp);
    phi_eta(1,pp) = .5+eta_lG(pp); 
    phi_eta(2,pp) = -.5+eta_lG(pp); 
    phi_eta(3,pp) = 0;
    phi_eta(4,pp) = 0; 
    phi_eta(5,pp) = -2*eta_lG(pp); 
    phi_eta(6,pp) = 0;
    for ii = 1:6
        for jj = 1:6
            Tmat2{pp}(ii,jj) = phi_xi(ii,pp)*phi_eta(jj,pp) ...
                              -phi_eta(ii,pp)*phi_xi(jj,pp);
        end
    end
end
save('Tmat2.mat','Tmat2')

%Creating matrix Tmat3
Tmat3 = cell(1,n_lGaussian_Q);
phi_xi = zeros(6,n_lGaussian_Q);
phi_eta = zeros(6,n_lGaussian_Q);
for pp = 1:n_lGaussian_Q
    phi_xi(1,pp) = 0; 
    phi_xi(2,pp) = .5; 
    phi_xi(3,pp) = .5+xi_lG(pp);
    phi_xi(4,pp) = 1-xi_lG(pp); 
    phi_xi(5,pp) = -1+xi_lG(pp); 
    phi_xi(6,pp) = -1-xi_lG(pp);
    phi_eta(1,pp) = .5-xi_lG(pp); 
    phi_eta(2,pp) = .5; 
    phi_eta(3,pp) = 0;
    phi_eta(4,pp) = 1+xi_lG(pp); 
    phi_eta(5,pp) = -1+xi_lG(pp); 
    phi_eta(6,pp) = -1-xi_lG(pp);
    for ii = 1:6
        for jj = 1:6
            Tmat3{pp}(ii,jj) = phi_xi(ii,pp)*phi_eta(jj,pp) ...
                              -phi_eta(ii,pp)*phi_xi(jj,pp);
        end
    end
end
save('Tmat3.mat','Tmat3')

%Creating matrix Tmat4
Tmat4 = cell(1,n_lGaussian_Q);
phi_xi = zeros(6,n_lGaussian_Q);
phi_eta = zeros(6,n_lGaussian_Q);
for pp = 1:n_lGaussian_Q
    phi_xi(1,pp) = 0; 
    phi_xi(2,pp) = .5; 
    phi_xi(3,pp) = .5+xi_lG(pp);
    phi_xi(4,pp) = 1-xi_lG(pp); 
    phi_xi(5,pp) = -1+xi_lG(pp); 
    phi_xi(6,pp) = -1-xi_lG(pp);
    phi_eta(1,pp) = .5-xi_lG(pp); 
    phi_eta(2,pp) = .5; 
    phi_eta(3,pp) = 0;
    phi_eta(4,pp) = 1+xi_lG(pp); 
    phi_eta(5,pp) = -1+xi_lG(pp); 
    phi_eta(6,pp) = -1-xi_lG(pp);
    for ii = 1:6
        for jj = 1:6
            Tmat4{pp}(ii,jj) = phi_xi(ii,pp)*phi_eta(jj,pp) ...
                              -phi_eta(ii,pp)*phi_xi(jj,pp);
        end
    end
end
save('Tmat4.mat','Tmat4')

%Creating matrix Tmat5
Tmat5 = cell(1,n_lGaussian_Q);
phi_xi = zeros(6,n_lGaussian_Q);
phi_eta = zeros(6,n_lGaussian_Q);
for pp = 1:n_lGaussian_Q
    phi_xi(1,pp) = 0; 
    phi_xi(2,pp) = -.5+eta_lG(pp); 
    phi_xi(3,pp) = -.5;
    phi_xi(4,pp) = 1+eta_lG(pp); 
    phi_xi(5,pp) = -1-eta_lG(pp); 
    phi_xi(6,pp) = 1-eta_lG(pp);
    phi_eta(1,pp) = .5+eta_lG(pp); 
    phi_eta(2,pp) = -.5+eta_lG(pp); 
    phi_eta(3,pp) = 0;
    phi_eta(4,pp) = 0; 
    phi_eta(5,pp) = -2*eta_lG(pp); 
    phi_eta(6,pp) = 0;
    for ii = 1:6
        for jj = 1:6
            Tmat5{pp}(ii,jj) = phi_xi(ii,pp)*phi_eta(jj,pp) ...
                              -phi_eta(ii,pp)*phi_xi(jj,pp);
        end
    end
end
save('Tmat5.mat','Tmat5')

%Line-element to triangular-element map for boundary 1
LE_to_TE_1 = zeros(n1_el,1);
for ee = 1:n1_el
    [R,~] = find(l1_1 == 2*ee);
    centre_node_global_num = l_1(R,2);
    [RT,~] = find(l == centre_node_global_num);
    LE_to_TE_1(ee) = RT;
end
save('LE_to_TE_1.mat','LE_to_TE_1')

%Line-element to triangular-element map for boundary 2
LE_to_TE_2 = zeros(n2_el,1);
for ee = 1:n2_el
    [R,~] = find(l2_2 == 2*ee);
    centre_node_global_num = l_2(R,2);
    [RT,~] = find(l == centre_node_global_num);
    LE_to_TE_2(ee) = RT;
end
save('LE_to_TE_2.mat','LE_to_TE_2')

%Line-element to triangular-element map for boundary 3
LE_to_TE_3 = zeros(n3_el,1);
for ee = 1:n3_el
    [R,~] = find(l3_3 == 2*ee);
    centre_node_global_num = l_3(R,2);
    [RT,~] = find(l == centre_node_global_num);
    LE_to_TE_3(ee) = RT;
end
save('LE_to_TE_3.mat','LE_to_TE_3')

%Line-element to triangular-element map for boundary 4
LE_to_TE_4 = zeros(n4_el,1);
for ee = 1:n4_el
    [R,~] = find(l4_4 == 2*ee);
    centre_node_global_num = l_4(R,2);
    [RT,~] = find(l == centre_node_global_num,1);
    LE_to_TE_4(ee) = RT;
end
save('LE_to_TE_4.mat','LE_to_TE_4')

%Line-element to triangular-element map for boundary 5
LE_to_TE_5 = zeros(n5_el,1);
for ee = 1:n5_el
    [R,~] = find(l5_5 == 2*ee);
    centre_node_global_num = l_5(R,2);
    [RT,~] = find(l == centre_node_global_num,1,'last');
    LE_to_TE_5(ee) = RT;
end
save('LE_to_TE_5.mat','LE_to_TE_5')

n_K_near = n_spines_near;
save('n_K_near.mat','n_K_near')
K_nodes_near = 1:n_spines_near;
save('K_nodes_near.mat','K_nodes_near')
n_K_far = n_spines_far;
save('n_K_far.mat','n_K_far')
K_nodes_far = n_spines_near+1:n_spines;
save('K_nodes_far.mat','K_nodes_far')

n_I_near = n_spines_near;
save('n_I_near.mat','n_I_near')
I_nodes_near = 1:n_spines_near;
save('I_nodes_near.mat','I_nodes_near')
n_I_far = n_spines_far;
save('n_I_far.mat','n_I_far')
I_nodes_far = n_spines_near+1:n_spines;
save('I_nodes_far.mat','I_nodes_far')

n_v_pre_start_far = n_v_start_far-1;
save('n_v_pre_start_far.mat','n_v_pre_start_far')

n_p_pre_start_far = n_p_start_far-1;
save('n_p_pre_start_far.mat','n_p_pre_start_far')

n_lambda2_uk = n_spines;
save('n_lambda2_uk.mat','n_lambda2_uk')

n_lambda2_uk_near = n_spines_near;
save('n_lambda2_uk_near.mat','n_lambda2_uk_near')
n_lambda2_uk_far = n_spines_far;
save('n_lambda2_uk_far.mat','n_lambda2_uk_far')

sep_nodes = n_v_start_far:n_v_near;
save('sep_nodes.mat','sep_nodes')

n_us2_uk_near = n_spines_near;
save('n_us2_uk_near.mat','n_us2_uk_near')
n_ws2_uk_near = n_spines_near;
save('n_ws2_uk_near.mat','n_ws2_uk_near')
n_rhos2_uk_near = n_spines_near;
save('n_rhos2_uk_near.mat','n_rhos2_uk_near')
n_sigma2_uk_near = n_spines_near;
save('n_sigma2_uk_near.mat','n_sigma2_uk_near')

n_us1_uk_near = n_spines_near;
save('n_us1_uk_near.mat','n_us1_uk_near')
n_ws1_uk_near = n_spines_near;
save('n_ws1_uk_near.mat','n_ws1_uk_near')
n_rhos1_uk_near = n_spines_near;
save('n_rhos1_uk_near.mat','n_rhos1_uk_near')
n_sigma1_uk_near = n_spines_near;
save('n_sigma1_uk_near.mat','n_sigma1_uk_near')

n_us2_uk_far = n_spines_far-1;
save('n_us2_uk_far.mat','n_us2_uk_far')
n_ws2_uk_far = n_spines_far;
save('n_ws2_uk_far.mat','n_ws2_uk_far')
n_rhos2_uk_far = n_spines_far;
save('n_rhos2_uk_far.mat','n_rhos2_uk_far')
n_sigma2_uk_far = n_spines_far;
save('n_sigma2_uk_far.mat','n_sigma2_uk_far')

n_us1_uk_far = n_spines_far-1;
save('n_us1_uk_far.mat','n_us1_uk_far')
n_ws1_uk_far = n_spines_far;
save('n_ws1_uk_far.mat','n_ws1_uk_far')
n_rhos1_uk_far = n_spines_far;
save('n_rhos1_uk_far.mat','n_rhos1_uk_far')
n_sigma1_uk_far = n_spines_far;
save('n_sigma1_uk_far.mat','n_sigma1_uk_far')


n_S2_near = n_spines_near;
save('n_S2_near.mat','n_S2_near')
n_E2_near = n_spines_near;
save('n_E2_near.mat','n_E2_near')
n_D2_near = n_spines_near;
save('n_D2_near.mat','n_D2_near')
n_T2_near = n_spines_near;
save('n_T2_near.mat','n_T2_near')

S2_nodes_near = 1:n_spines_near;
save('S2_nodes_near.mat','S2_nodes_near')
E2_nodes_near = 1:n_spines_near;
save('E2_nodes_near.mat','E2_nodes_near')
D2_nodes_near = 1:n_spines_near;
save('D2_nodes_near.mat','D2_nodes_near')
T2_nodes_near = 1:n_spines_near;
save('T2_nodes_near.mat','T2_nodes_near')

n_S1_near = n_spines_near;
save('n_S1_near.mat','n_S1_near')
n_E1_near = n_spines_near;
save('n_E1_near.mat','n_E1_near')
n_D1_near = n_spines_near;
save('n_D1_near.mat','n_D1_near')
n_T1_near = n_spines_near;
save('n_T1_near.mat','n_T1_near')

S1_nodes_near = 1:n_spines_near;
save('S1_nodes_near.mat','S1_nodes_near')
E1_nodes_near = 1:n_spines_near;
save('E1_nodes_near.mat','E1_nodes_near')
D1_nodes_near = 1:n_spines_near;
save('D1_nodes_near.mat','D1_nodes_near')
T1_nodes_near = 1:n_spines_near;
save('T1_nodes_near.mat','T1_nodes_near')

n_S2_far = n_spines_far;
save('n_S2_far.mat','n_S2_far')
n_E2_far = n_spines_far;
save('n_E2_far.mat','n_E2_far')
n_D2_far = n_spines_far-1;
save('n_D2_far.mat','n_D2_far')
n_T2_far = n_spines_far;
save('n_T2_far.mat','n_T2_far')

S2_nodes_far = n_spines_near+1:n_spines;
save('S2_nodes_far.mat','S2_nodes_far')
E2_nodes_far = n_spines_near+1:n_spines;
save('E2_nodes_far.mat','E2_nodes_far')
D2_nodes_far = n_spines_near+1:n_spines-1;
save('D2_nodes_far.mat','D2_nodes_far')
T2_nodes_far = n_spines_near+1:n_spines;
save('T2_nodes_far.mat','T2_nodes_far')

n_S1_far = n_spines_far;
save('n_S1_far.mat','n_S1_far')
n_E1_far = n_spines_far;
save('n_E1_far.mat','n_E1_far')
n_D1_far = n_spines_far-1;
save('n_D1_far.mat','n_D1_far')
n_T1_far = n_spines_far;
save('n_T1_far.mat','n_T1_far')

S1_nodes_far = n_spines_near+1:n_spines;
save('S1_nodes_far.mat','S1_nodes_far')
E1_nodes_far = n_spines_near+1:n_spines;
save('E1_nodes_far.mat','E1_nodes_far')
D1_nodes_far = n_spines_near+1:n_spines-1;
save('D1_nodes_far.mat','D1_nodes_far')
T1_nodes_far = n_spines_near+1:n_spines;
save('T1_nodes_far.mat','T1_nodes_far')

lambda2nodes_nums_uk_near = 1:n_spines_near;
save('lambda2nodes_nums_uk_near.mat','lambda2nodes_nums_uk_near')
us2nodes_nums_uk_near = 1:n_spines_near;
save('us2nodes_nums_uk_near.mat','us2nodes_nums_uk_near')
ws2nodes_nums_uk_near = 1:n_spines_near;
save('ws2nodes_nums_uk_near.mat','ws2nodes_nums_uk_near')
rhos2nodes_nums_uk_near = 1:n_spines_near;
save('rhos2nodes_nums_uk_near.mat','rhos2nodes_nums_uk_near')
sigma2nodes_nums_uk_near = 1:n_spines_near;
save('sigma2nodes_nums_uk_near.mat','sigma2nodes_nums_uk_near')

spine_nums_uk_near = 1:n_spines_near;
save('spine_nums_uk_near.mat','spine_nums_uk_near')
us1nodes_nums_uk_near = 1:n_spines_near;
save('us1nodes_nums_uk_near.mat','us1nodes_nums_uk_near')
ws1nodes_nums_uk_near = 1:n_spines_near;
save('ws1nodes_nums_uk_near.mat','ws1nodes_nums_uk_near')
rhos1nodes_nums_uk_near = 1:n_spines_near;
save('rhos1nodes_nums_uk_near.mat','rhos1nodes_nums_uk_near')
sigma1nodes_nums_uk_near = 1:n_spines_near;
save('sigma1nodes_nums_uk_near.mat','sigma1nodes_nums_uk_near')

lambda2nodes_nums_uk_far = n_spines_near+1:n_spines;
save('lambda2nodes_nums_uk_far.mat','lambda2nodes_nums_uk_far')
us2nodes_nums_uk_far = n_spines_near+1:n_spines-1;
save('us2nodes_nums_uk_far.mat','us2nodes_nums_uk_far')
ws2nodes_nums_uk_far = n_spines_near+1:n_spines;
save('ws2nodes_nums_uk_far.mat','ws2nodes_nums_uk_far')
rhos2nodes_nums_uk_far = n_spines_near+1:n_spines;
save('rhos2nodes_nums_uk_far.mat','rhos2nodes_nums_uk_far')
sigma2nodes_nums_uk_far = n_spines_near+1:n_spines;
save('sigma2nodes_nums_uk_far.mat','sigma2nodes_nums_uk_far')

spine_nums_uk_far = n_spines_near+1:n_spines;
save('spine_nums_uk_far.mat','spine_nums_uk_far')
us1nodes_nums_uk_far = n_spines_near+1:n_spines-1;
save('us1nodes_nums_uk_far.mat','us1nodes_nums_uk_far')
ws1nodes_nums_uk_far = n_spines_near+1:n_spines;
save('ws1nodes_nums_uk_far.mat','ws1nodes_nums_uk_far')
rhos1nodes_nums_uk_far = n_spines_near+1:n_spines;
save('rhos1nodes_nums_uk_far.mat','rhos1nodes_nums_uk_far')
sigma1nodes_nums_uk_far = n_spines_near+1:n_spines;
save('sigma1nodes_nums_uk_far.mat','sigma1nodes_nums_uk_far')

c1 = 1; %Node number of the contact line in the boundary 1
%local numbering system
save('c1.mat','c1')
c2 = 1; %Node number of the contact line in the boundary 2
%local numbering system
save('c2.mat','c2')

a1 = n_spines; %Node number of the apex in the boundary 1 numbering system
save('a1.mat','a1')

o2 = n_spines; %Node number of the origin in the boundary 2 numbering system
save('o2.mat','o2')

press_lim_node_solid = 2;
save('press_lim_node_solid.mat','press_lim_node_solid')
press_lim_node_gas = 2 + n_el_sing;
save('press_lim_node_gas.mat','press_lim_node_gas')

n_spines_pre_start_far = n_spines_near-1;
save('n_spines_pre_start_far.mat','n_spines_pre_start_far')

%number of the far-field boundary-1-node residual where the joint is
%located
%this is the number of the far field equation that has to be moved to the
%near field
b1Jeq = 1;
save('b1Jeq.mat','b1Jeq')

%number of the far-field boundary-2-node residual where the joint is
%located
%this is the number of the far field equation that has to be moved to the
%near field
b2Jeq = 1;
save('b2Jeq.mat','b2Jeq')
