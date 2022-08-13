close all
clear
clc
 
% Dimensional fluid constants (tin as leaves cartridge) for NS
g_dim = 9.8; %m/s^2
save('g_dim.mat','g_dim') 
rho_dim = 1000; %kg/m^3
save('rho_dim.mat','rho_dim') 
nu_dim = 1E-6; %m^2/s
save('nu_dim.mat','nu_dim') 
sigma_1e_dim = 72.86E-3;%Kg/s^2 (i.e. N/m)
save('sigma_1e_dim.mat','sigma_1e_dim')
R = .5E-6; %droplet radius

%equilibrium contact angle
theta_c_eq = 60*pi/180;
save('theta_c_eq.mat','theta_c_eq')

% width of the interfacial phase in metres (recommended 
%value is 1E-9)
width_dim = 1E-9; 
save('width_dim.mat','width_dim')

%Data for initial configuration of the simulation 
theta_c_ini = 160*pi/180; %radians 
save('theta_c_ini.mat','theta_c_ini')

%factor by which we devide the minimum lengthof the mesh that was
%recommended by James
factor = 2^12;

%spine feet parameters
%this part here sets the number and spacing of the spines closer to the
%contact line. In this region, each odd-numbered spine has one more 
%element to its right than to its left
spine_factor_r =  1.01;%1.15;%;%needs to be slightly bigger than 1     1.05;%
save('spine_factor_r.mat','spine_factor_r')
n_spines_r = 5;% 5;% must be an odd number          5;%
save('n_spines_r.mat','n_spines_r')

spine_factor_mr = 1.01;%                   2;%   
save('spine_factor_mr.mat','spine_factor_mr')
n_spines_mr = 3;% 5;% must be an odd number        
save('n_spines_mr.mat','n_spines_mr')

%This part defines the number of spines and their separation for the
%following section of the domain. Spines can be more separated here as we 
%do not need such a high resultion furhter away from the contact line
spine_factor_m = 1.01;%                   2;%   
save('spine_factor_m.mat','spine_factor_m')
n_spines_m = 3;% 5;% must be an odd number        
save('n_spines_m.mat','n_spines_m') 

%The following part sets the number and spacing of the spines near the
%inflow boundary. An intermediate section between this one and the last one
%constructed is then created using uniform spacing
l_min_l = 5E-3 ;% 1E-2;%               .01;%
spine_factor_l = 1.2;% 1.11;% this needs to be slightly bigger than 1, James' value = 1.07
save('spine_factor_l.mat','spine_factor_l')
n_spines_l = 3;% 5;% must be an odd number        5;% 
save('n_spines_l.mat','n_spines_l')

%number of singular elements at the contact line
n_el_sing = 3;
save('n_el_sing.mat','n_el_sing')

%% IFM coefficients for the free-surface
%Beta for liquid-gas interface (this
%could have a coefficient in front, but we take it to be one)
beta_g_dim = rho_dim*nu_dim/width_dim;
save('beta_g_dim.mat','beta_g_dim')
%Alpha coefficient for liquid-gas interface (this could also
%have a coefficient in front but we also take it to be one, 
%making it the inverse of the beta_g constant)
alpha_g_dim = width_dim/(rho_dim*nu_dim); 
save('alpha_g_dim.mat','alpha_g_dim')
%Surface density at which surface tension is zero 
%(this is taken as the volume density times the width of the 
%surface phase)
rhos_0_dim = rho_dim*width_dim; 
save('rhos_0_dim.mat','rhos_0_dim')
%Equilibrium surface density of liquid-gas 
%interface
rhos1_e_dim = 0.6*rhos_0_dim; 
save('rhos1_e_dim.mat','rhos1_e_dim')
%Proportionality constant for the state equation for surface 
%tension as a function of density on the liquid-gas 
%interface
gamma_g_dim = sigma_1e_dim/(rhos_0_dim-rhos1_e_dim);
save('gamma_g_dim.mat','gamma_g_dim')
%relaxation time of surface 1 in seconds
tau_g_dim = rho_dim*nu_dim*7E-6;
save('tau_g_dim.mat','tau_g_dim')

%solid-gas surface tension
sigma_gs = 0;

%IFM coefficients for the liquid-solid interface
%constant in front of the state equation for the liquid-solid
%surface
gamma_s_dim = gamma_g_dim;
save('gamma_s_dim.mat','gamma_s_dim')
%Equilibrium surface tension for the liqud-solid interface
%(this is a function of other paremeter that we have already
%chosen)
sigma_2e_dim = -sigma_1e_dim*cos(theta_c_eq)+sigma_gs;
save('sigma_2e_dim.mat','sigma_2e_dim')
%equilibrium surface density for liquid-solid surface
%(chosen so as to have a 30 degree static contact angle)
rhos2_e_dim = rhos_0_dim - sigma_2e_dim/gamma_s_dim; 
save('rhos2_e_dim.mat','rhos2_e_dim')
%Beta navier slip
beta_ns_dim = rho_dim*nu_dim/width_dim;   
save('beta_ns_dim.mat','beta_ns_dim')
%alpha coefficient for the solid-liquid interface
alpha_s_dim = alpha_g_dim; 
save('alpha_s_dim.mat','alpha_s_dim')
%relaxation time for the liquid-solid interface
tau_s_dim = rho_dim*nu_dim*7E-6;
save('tau_s_dim.mat','tau_s_dim')

%Unit height 
L_dim = R
save('L_dim.mat','L_dim') 

%Unit speed
U_dim = sigma_1e_dim/(rho_dim*nu_dim) %i.e. Ca = 1
% U_dim = .01 %m/s
save('U_dim.mat','U_dim')

%Unit time (capillary time based on unit length)
T_dim = L_dim/U_dim 
save('T_dim.mat','T_dim')

%Unit of pressure 
P_dim = rho_dim*nu_dim*U_dim/L_dim;
save('P_dim.mat','P_dim')

%Dimensionless numbers for the bulk equations
Re =  U_dim*L_dim/nu_dim %Reynolds number
save('Re.mat','Re')
St = 0 %g_dim*L_dim^2/(nu_dim*U_dim); %Stokes number
save('St.mat','St')
Ca = rho_dim*nu_dim*U_dim/sigma_1e_dim %Capillary number
save('Ca.mat','Ca')

%Dimensionless numbers for the liquid-gas IFM
%dimensionless beta coefficient for the gas interface
Bg = beta_g_dim*U_dim*L_dim/sigma_1e_dim
save('Bg.mat','Bg')
%dimensionless alpha coefficient for the gas interface
Eg = alpha_g_dim*sigma_1e_dim/(U_dim*L_dim)
save('Eg.mat','Eg')
%Dimensionless constant for the state equation on
%the liquid-gas interface
Cg = gamma_g_dim*rhos_0_dim/sigma_1e_dim
save('Cg.mat','Cg')
%dimensionless equilibrium density
Dg = rhos1_e_dim/rhos_0_dim
save('Dg.mat','Dg')
%dimensionless surface density ratio
Lg = width_dim/L_dim
save('Lg.mat','Lg')
% Fg = rhos_0_dim/(rho_dim*U_dim*tau_g_dim)
% save('Fg.mat','Fg')
%Dimensionless relaxation time
Tg = tau_g_dim*U_dim/L_dim
save('Tg.mat','Tg')

%Dimensionless numbers for liquid-gas IFM
%dimensionless beta coefficient for the gas interface
Be = beta_ns_dim*L_dim/(rho_dim*nu_dim)
save('Be.mat','Be')
%dimensionless alpha coefficient for the solid-liquid 
%interface
Es = alpha_s_dim*sigma_1e_dim/(U_dim*L_dim)
save('Es.mat','Es')
%dimensionless coefficient of the state equation for the
%liquid-solid interface
Cs = gamma_s_dim*rhos_0_dim/sigma_1e_dim
save('Cs.mat','Cs')
%dimensionless equilibrium density
Ds = rhos2_e_dim/rhos_0_dim
save('Ds.mat','Ds')
%dimensionless surface-density ratio for the liquid-solid
%interface
Ls = width_dim/L_dim;
save('Ls.mat','Ls')
% Fs = rhos_0_dim/(rho_dim*U_dim*tau_s_dim)
% save('Fs.mat','Fs')
%Dimensionless relaxation time for solid-liquid interface
Ts = tau_s_dim*U_dim/L_dim
save('Ts.mat','Ts')

%Contact line dimensionless numbers
So = sigma_gs/sigma_1e_dim
save('So.mat','So')

%z-traslation for the centre of the sphere
trans_z = cos(pi-theta_c_ini);

%distance from origin to apex
h_0 = 1+trans_z;

%initial distance from origin to focus
f = sin(pi-theta_c_ini); %based on the shape of a parabola

%Required mesh resolution
l_min = 0.03;% min([5/(100*Bg) 1/Be  Ts*Ca/10])/(factor );% .05;% following J&Y IJFNMF 2012 

tol_spine_ang = 1E-12; %this is just a numerical parameter to define 
%convergence of the spine lentgh's initial guess

%% Creating spine feet

spine_feet = zeros(1,n_spines_r+n_spines_m-1);
spine_feet(1) = f;
spine_feet(2) = f-l_min;
spine_feet(3) = f-l_min*(1+spine_factor_r);
for ii = 4:n_spines_r
    spine_feet(ii) = spine_feet(ii-1) ...
                     - spine_factor_r*(spine_feet(ii-2) ...
                                       - spine_feet(ii-1));
end
figure
plot(spine_feet(1:n_spines_r),zeros(1,n_spines_r),'xb')
hold on
for ii = n_spines_r+1:n_spines_r+n_spines_mr-1
    spine_feet(ii) = spine_feet(ii-1) ...
                     - spine_factor_mr*(spine_feet(ii-2) ...
                                       - spine_feet(ii-1));
end
plot(spine_feet(n_spines_r+1:n_spines_r+n_spines_mr-1),zeros(1,n_spines_mr-1),'^k')
for ii = n_spines_r+n_spines_mr:n_spines_r+n_spines_mr+n_spines_m-2
    spine_feet(ii) = spine_feet(ii-1) ...
                     - spine_factor_m*(spine_feet(ii-2) ...
                                       - spine_feet(ii-1));
end
plot(spine_feet(n_spines_r+n_spines_mr:n_spines_r+n_spines_mr+n_spines_m-2),zeros(1,n_spines_m-1),'+g')
plot(0,0,'db')
pause
%constructing the spine feet on the left end
spines_feet_left_end = zeros(1,n_spines_l);
spines_feet_left_end(end) = 0;
spines_feet_left_end(end-1) = 0 + l_min_l;
for ii = 2:n_spines_l-1
    spines_feet_left_end(end-ii) = spines_feet_left_end(end-ii+1) ...
                                   + spine_factor_l ...
                                     *(spines_feet_left_end(end-ii+1) ...
                                       - spines_feet_left_end(end-ii+2));
end
plot(spines_feet_left_end,zeros(1,n_spines_l),'xg')
pause
left_width = spines_feet_left_end(1)-spines_feet_left_end(2);
right_width = spine_feet(end-1)-spine_feet(end);
av_width = (left_width+right_width)/2;
n_inter = ceil((spine_feet(end)-spines_feet_left_end(1))/av_width);
if mod(n_inter,2) ~= 0
    n_inter = n_inter+1;
end
width = (spine_feet(end)-spines_feet_left_end(1))/n_inter;
spine_feet = [spine_feet,zeros(1,n_inter-1)];
for ii = 0:n_inter-2
    spine_feet(end-ii) = spines_feet_left_end(1)+(ii+1)*width;
end
plot(spine_feet(end-n_inter+2:end),zeros(1,n_inter-1),'dk')
pause
%Joining everything
spine_feet = [spine_feet,spines_feet_left_end];
n_spines = length(spine_feet);
save('n_spines.mat','n_spines')
close all
figure
plot(spine_feet,zeros(size(spine_feet)),'xb');

save('spine_feet.mat','spine_feet')
spine_feet_locator = spine_feet/f;
save('spine_feet_locator.mat','spine_feet_locator')

%Chi value for each spine
%This is the coordinate of the bipolar system that corresponds to the
%circle that passes through each spine feet
chi_spines = zeros(1,n_spines);
chi_spines(1) = Inf;
chi_spines(end) = 0;
for ii = 2:n_spines-1
    fun = @(chi) foot_of_chi(chi,f)-spine_feet(ii);
    chi_spines(ii) = fzero(fun,5);
end

%The radius of each spine (from the third spine on) is determined by its
%chi coordinate
R_spines = zeros(1,n_spines);
%we force the second spine to have its centre at the contact line. This is
%so that we can enforce the condition lim_{r\to 0} \partial_{\theta}p = 0,
%by equating the pressure values on the two points of this spine as done in
%S&S CMAME 2011
R_spines(2) = abs(spine_feet(1)-spine_feet(2));
R_spines(3) = abs(spine_feet(1)-spine_feet(3));
R_spines(end) = inf;
for ii = 4:n_spines-1
    R_spines(ii) = f/sinh(chi_spines(ii));
end
save('R_spines.mat','R_spines')
 
%spine length vector intial form
spine_lengths_ini = zeros(n_spines,1);
spine_lengths_ini(1) = f; %this spine is always zero, so we just store f here
spine_lengths_ini(end) = h_0;
%Here we determine the length of each spine (so as to match the initial
%guess for our surface profile)
for ii = 2:n_spines-1
    %We use the bisection method on the angle that each spine defines
    %with a tolerance pre-defined above to assess whet
    theta_l = 0;
%     if ii < 80
        theta_u = 1.01*theta_c_ini;
%     else
%         theta_u = asin(2/R_spines(ii));
%     end
    x_l = spine_feet(ii) + (1-cos(theta_l))*R_spines(ii);
    y_l = sin(theta_l)*R_spines(ii);
    %function f is positive outside the curve that defines the initial
    %surface profile (a circle) and negative inside it 
    f_l = x_l^2+(y_l-cos(pi-theta_c_ini))^2-1;
    x_u = spine_feet(ii) + (1-cos(theta_u))*R_spines(ii);
    y_u = sin(theta_u)*R_spines(ii);
    f_u = x_u^2+(y_u-cos(pi-theta_c_ini))^2-1;
    if f_l*f_u > 0
         disp('warning: wrong initialisation of bisection method') 
         pause
    end
    while theta_u-theta_l > tol_spine_ang
        theta_m = (theta_l + theta_u)/2;
        x_m = spine_feet(ii) + (1-cos(theta_m))*R_spines(ii);
        y_m = sin(theta_m)*R_spines(ii);
        f_m = x_m^2+(y_m-cos(pi-theta_c_ini))^2-1;
        if f_l*f_m > 0
            theta_l = theta_m;
            f_l = f_m;
        elseif f_l*f_m < 0
            theta_u = theta_m;
            f_u = f_m;
        else
            f_l = f_m;
            f_u = f_m;
            theta_l = theta_m;
            theta_u = theta_m;
        end
    end
    spine_lengths_ini(ii) = theta_m*R_spines(ii);
end
save('spine_lengths_ini.mat','spine_lengths_ini')

%% Creating nodes

%%Finding total number of global nodes 

%First, we need to decide on which spines we will have an increasing number
%of nodes and on which we will have a fixed number of nodes
%The code is set-up so that the first spines have an increasing number of
%nodes. More specifically, from the 4th spine on (including the 4th spine)
%each spine has one more velocity node that the previous spine, and this
%behaviour goes on until spine number n_spines_incr, defined below. The
%code is set up so that n_spines_incr can be any odd number between
%n_spines_r and n_spines_r+n_spines_mr-1
n_spines_incr = n_spines_r;%+n_spines_mr-1;
save('n_spines_incr.mat','n_spines_incr')

n_spines_const = n_spines-n_spines_incr;
save('n_spines_const.mat','n_spines_const')

n_nodes_const_per_spine = 2*n_el_sing + n_spines_incr - 2; %number of nodes 
%in each of the spines with constant nodes
save('n_nodes_const_per_spine.mat','n_nodes_const_per_spine')

%Sum of the number in the first spine (1), plus an arithmetic progression 
%plus a constant times the number of spines with a fixed number of nodes on
%them
n_v_incr = 2 + n_el_sing ...
           + (2*n_el_sing + 1 + 2*n_el_sing + n_spines_incr - 2) ...
             *(n_spines_incr-2)/2;
save('n_v_incr.mat','n_v_incr')

n_v = n_v_incr ...
      + n_nodes_const_per_spine*n_spines_const;
save('n_v.mat','n_v')

%nodes (i.e. matrix containing the (r,z) coordinates of all nodes where 
%each column is a node and they are listed in order)
Nodes_rz = zeros(n_v,2);
Nodes_rz(1,:) = [f ,0]; 
alpha_s = zeros(n_v,1); %This vector contains the data on how far along its
%spine the node is
Nodes_rz(2,1) = spine_feet(2);
Nodes_rz(2,2) = 0;
for kk = 3:n_el_sing+2
    alpha_s(kk) = (kk-2)/n_el_sing;
    %angle in the constant chi circle for each node 
    %(arc length over radius)
    theta = alpha_s(kk)*spine_lengths_ini(2)/R_spines(2);
    %x of the node is the foot of the spine plus the circle's radius 
    %minus (radius*cos(theta))
    Nodes_rz(kk,1) = spine_feet(2) + R_spines(2)*(1-cos(theta));
    Nodes_rz(kk,2) = R_spines(2)*sin(theta);
end
for ii = 3:n_spines_incr
    kk = kk+1;
    Nodes_rz(kk,1) = spine_feet(ii);
    Nodes_rz(kk,2) = 0;
    for jj = 2:2*n_el_sing+ii-2
        kk = kk+1;
        alpha_s(kk) = (jj-1)/(2*n_el_sing+ii-3);
        %angle in the constant chi circle for each node 
        %(arc length over radius)
        theta = alpha_s(kk)*spine_lengths_ini(ii)/R_spines(ii); 
        %x of the node is the foot of the spine plus the circle's radius 
        %minus (radius*cos(theta))
        Nodes_rz(kk,1) = spine_feet(ii) + R_spines(ii)*(1-cos(theta));
        Nodes_rz(kk,2) = R_spines(ii)*sin(theta);
    end
end
for ii = n_spines_incr+1:n_spines-1
    kk = kk+1;
    Nodes_rz(kk,1) = spine_feet(ii);
    Nodes_rz(kk,2) = 0;
    for jj = 2:n_nodes_const_per_spine
        kk = kk+1;
        alpha_s(kk) = (jj-1)/(n_nodes_const_per_spine-1);
        %angle in the constant chi circle for each node 
        %(arc length over radius)
        theta = alpha_s(kk)*spine_lengths_ini(ii)/R_spines(ii); 
        %x of the node is the foot of the spine plus the circle's radius 
        %minus the radius cos(theta)
        Nodes_rz(kk,1) = spine_feet(ii) + R_spines(ii)*(1-cos(theta));
        Nodes_rz(kk,2) = R_spines(ii)*sin(theta);
    end
end
for ii = n_v-n_nodes_const_per_spine+1:n_v 
    alpha_s(ii) = (ii-(n_v-n_nodes_const_per_spine+1)) ...
                  /(n_nodes_const_per_spine-1);
end 
save('alpha_s.mat','alpha_s')
Nodes_rz(end-n_nodes_const_per_spine+1:end,1) = zeros(1,n_nodes_const_per_spine);
Nodes_rz(end-n_nodes_const_per_spine+1:end,2) = 0:h_0/(n_nodes_const_per_spine-1):h_0;
save('Nodes_rz.mat','Nodes_rz')
 


n_v %#ok
pause
close all
%Nodes
figure
plot(Nodes_rz(:,1),Nodes_rz(:,2),'x')
axis equal

%Spine feet and tip node numbers
spine_feet_nodes = ones(1,n_spines);
spine_tip_nodes = ones(1,n_spines);
spine_feet_nodes(2) = 2;
spine_tip_nodes(2) = 2 + n_el_sing;
for ii = 3:n_spines_incr
    spine_feet_nodes(ii) = spine_tip_nodes(ii-1)+1;
    spine_tip_nodes(ii) = 2 + n_el_sing ...
                          + (2*n_el_sing + 1 + 2*n_el_sing + ii - 2) ...
                            *(ii-2)/2;
end
%the spines in the region where spines have a fixed number of nodes.
for ii = n_spines_incr+1:n_spines
    spine_feet_nodes(ii) = spine_tip_nodes(ii-1)+1;
    spine_tip_nodes(ii) = spine_tip_nodes(ii-1) + n_nodes_const_per_spine;
end
save('spine_feet_nodes.mat','spine_feet_nodes')
save('spine_tip_nodes.mat','spine_tip_nodes')



