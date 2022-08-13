close all
clear
clc

load('Nodes_rz.mat')
load('l.mat')
load('n_el.mat')
load('l_1.mat')
load('n1_el.mat')
load('l_2.mat')
load('n2_el.mat')
load('l_3.mat')
load('n3_el.mat')
load('n_v.mat');
load('lp.mat');
load('l1_1.mat')
load('l2_2.mat')
load('l3_3.mat')
load('n_spines.mat')
load('n_spines_incr.mat')
load('spine_feet.mat')
load('spine_feet_locator.mat')
load('spine_lengths_ini.mat')
load('spine_feet_nodes.mat')
load('spine_tip_nodes.mat')
load('n_el_near.mat')
load('alpha_s.mat')
load('n_nodes_const_per_spine.mat')
load('n_el_sing.mat')
f = spine_lengths_ini(1);
n_smooth = 128;
n_smooth_spines = 256;

% vector used for the curved sides
tvec = -1:2/n_smooth:1;
phi_l = (tvec-1).*tvec/2;
phi_c = (tvec+1).*(1-tvec);
phi_r = (tvec+1).*tvec/2;








%Domain boundaries
h0 = figure;
hold on
for ii = 1:n1_el
    xvec = Nodes_rz(l_1(ii,1),1)*phi_l+Nodes_rz(l_1(ii,2),1)*phi_c+Nodes_rz(l_1(ii,3),1)*phi_r;
    yvec = Nodes_rz(l_1(ii,1),2)*phi_l+Nodes_rz(l_1(ii,2),2)*phi_c+Nodes_rz(l_1(ii,3),2)*phi_r;
    plot(xvec,yvec,'b','LineWidth',2)
end
for ii = 1:n2_el
    xvec = Nodes_rz(l_2(ii,1),1)*phi_l+Nodes_rz(l_2(ii,2),1)*phi_c+Nodes_rz(l_2(ii,3),1)*phi_r;
    yvec = Nodes_rz(l_2(ii,1),2)*phi_l+Nodes_rz(l_2(ii,2),2)*phi_c+Nodes_rz(l_2(ii,3),2)*phi_r;
    plot(xvec,yvec,'g','LineWidth',2)
end
for ii = 1:n3_el
    xvec = Nodes_rz(l_3(ii,1),1)*phi_l+Nodes_rz(l_3(ii,2),1)*phi_c+Nodes_rz(l_3(ii,3),1)*phi_r;
    yvec = Nodes_rz(l_3(ii,1),2)*phi_l+Nodes_rz(l_3(ii,2),2)*phi_c+Nodes_rz(l_3(ii,3),2)*phi_r;
    plot(xvec,yvec,'r','LineWidth',2)
end
title('Domain and boundaries')
set(gca,'Fontsize',16)
text(.35,1,'Liquid','Fontsize',16)
text(1.5,.8,'Gas','Fontsize',16)
axis equal
set(gca,'xlim',[-.01 2],'ylim',[-.01 2.01])
box on
print('-depsc','Domain.eps')







%Domain, spines and nodes
h1 = figure;
hold on
%plotting boundary 1
for ii = 1:n1_el
    xvec = Nodes_rz(l_1(ii,1),1)*phi_l+Nodes_rz(l_1(ii,2),1)*phi_c+Nodes_rz(l_1(ii,3),1)*phi_r;
    yvec = Nodes_rz(l_1(ii,1),2)*phi_l+Nodes_rz(l_1(ii,2),2)*phi_c+Nodes_rz(l_1(ii,3),2)*phi_r;
    plot(xvec,yvec,'b','LineWidth',2)
end
%plotting boundary 2
for ii = 1:n2_el
    xvec = Nodes_rz(l_2(ii,1),1)*phi_l+Nodes_rz(l_2(ii,2),1)*phi_c+Nodes_rz(l_2(ii,3),1)*phi_r;
    yvec = Nodes_rz(l_2(ii,1),2)*phi_l+Nodes_rz(l_2(ii,2),2)*phi_c+Nodes_rz(l_2(ii,3),2)*phi_r;
    plot(xvec,yvec,'g','LineWidth',2)
end
%plotting boundary 3
for ii = 1:n3_el
    xvec = Nodes_rz(l_3(ii,1),1)*phi_l+Nodes_rz(l_3(ii,2),1)*phi_c+Nodes_rz(l_3(ii,3),1)*phi_r;
    yvec = Nodes_rz(l_3(ii,1),2)*phi_l+Nodes_rz(l_3(ii,2),2)*phi_c+Nodes_rz(l_3(ii,3),2)*phi_r;
    plot(xvec,yvec,'r','LineWidth',2)
end
f = spine_lengths_ini(1);
%Chi value for each spine
chi_spines = zeros(1,n_spines);
chi_spines(1) = Inf;
chi_spines(end) = 0;
for ii = 2:n_spines-1
    fun = @(chi) foot_of_chi(chi,f)-spine_feet(ii);
    chi_spines(ii) = fzero(fun,5);
end
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
%plotting spines
theta = zeros(1,n_smooth_spines);
xs = zeros(1,n_smooth_spines);
ys = zeros(1,n_smooth_spines);
for ss = 2:n_spines-1
    for jj = 1:n_smooth_spines
        %angle in the constant chi circle for each node (arc length over radius)
        theta(jj) = ((jj-1)/(n_smooth_spines-1))*spine_lengths_ini(ss)/R_spines(ss); 
        %x of the node is the foot of the spine plus the circle's radius 
        %minus the radius cos(theta)
        xs(jj) = spine_feet(ss) + R_spines(ss)*(1-cos(theta(jj)));
        ys(jj) = R_spines(ss)*sin(theta(jj));
    end
    plot(xs,ys,'-k')
end
%plotting the velocity nodes
for ii = 1:n_el
    plot(Nodes_rz(l(ii,:),1),Nodes_rz(l(ii,:),2),'ok','MarkerSize',6,'MarkerFaceColor','k')
end
%plotting global node numbers
for ii = 1:n_v
    text(Nodes_rz(ii,1),Nodes_rz(ii,2),[' i = ',num2str(ii)],'FontSize',24)
end
title('Spines and nodes')
set(gca,'Fontsize',24)










% Mesh
h2_a = figure;
hold on
for ii = n_el_near+1:n_el
    xvec = Nodes_rz(l(ii,1),1)*phi_l+Nodes_rz(l(ii,5),1)*phi_c+Nodes_rz(l(ii,2),1)*phi_r;
    yvec = Nodes_rz(l(ii,1),2)*phi_l+Nodes_rz(l(ii,5),2)*phi_c+Nodes_rz(l(ii,2),2)*phi_r;
    plot(xvec,yvec,'-k')
    xvec = Nodes_rz(l(ii,2),1)*phi_l+Nodes_rz(l(ii,6),1)*phi_c+Nodes_rz(l(ii,3),1)*phi_r;
    yvec = Nodes_rz(l(ii,2),2)*phi_l+Nodes_rz(l(ii,6),2)*phi_c+Nodes_rz(l(ii,3),2)*phi_r;
    plot(xvec,yvec,'-k')
    xvec = Nodes_rz(l(ii,3),1)*phi_l+Nodes_rz(l(ii,4),1)*phi_c+Nodes_rz(l(ii,1),1)*phi_r;
    yvec = Nodes_rz(l(ii,3),2)*phi_l+Nodes_rz(l(ii,4),2)*phi_c+Nodes_rz(l(ii,1),2)*phi_r;
    plot(xvec,yvec,'-k')
end
set(gca,'FontSize',24)
axis equal
title('Mesh')

% Mesh
h2_b = figure;
hold on
for ii = 1:n_el_near
    xvec = Nodes_rz(l(ii,1),1)*phi_l+Nodes_rz(l(ii,5),1)*phi_c+Nodes_rz(l(ii,2),1)*phi_r;
    yvec = Nodes_rz(l(ii,1),2)*phi_l+Nodes_rz(l(ii,5),2)*phi_c+Nodes_rz(l(ii,2),2)*phi_r;
    plot(xvec,yvec,'-k')
    xvec = Nodes_rz(l(ii,2),1)*phi_l+Nodes_rz(l(ii,6),1)*phi_c+Nodes_rz(l(ii,3),1)*phi_r;
    yvec = Nodes_rz(l(ii,2),2)*phi_l+Nodes_rz(l(ii,6),2)*phi_c+Nodes_rz(l(ii,3),2)*phi_r;
    plot(xvec,yvec,'-k')
    xvec = Nodes_rz(l(ii,3),1)*phi_l+Nodes_rz(l(ii,4),1)*phi_c+Nodes_rz(l(ii,1),1)*phi_r;
    yvec = Nodes_rz(l(ii,3),2)*phi_l+Nodes_rz(l(ii,4),2)*phi_c+Nodes_rz(l(ii,1),2)*phi_r;
    plot(xvec,yvec,'-k')
end
set(gca,'FontSize',24)
axis equal
title('Mesh')






% 
% 
% 
% 
%mesh with element numbers
h2 = figure;
hold on
for ii = 1:n_el
    xvec = Nodes_rz(l(ii,1),1)*phi_l+Nodes_rz(l(ii,5),1)*phi_c+Nodes_rz(l(ii,2),1)*phi_r;
    yvec = Nodes_rz(l(ii,1),2)*phi_l+Nodes_rz(l(ii,5),2)*phi_c+Nodes_rz(l(ii,2),2)*phi_r;
    plot(xvec,yvec,'-k')
    xvec = Nodes_rz(l(ii,2),1)*phi_l+Nodes_rz(l(ii,6),1)*phi_c+Nodes_rz(l(ii,3),1)*phi_r;
    yvec = Nodes_rz(l(ii,2),2)*phi_l+Nodes_rz(l(ii,6),2)*phi_c+Nodes_rz(l(ii,3),2)*phi_r;
    plot(xvec,yvec,'-k')
    xvec = Nodes_rz(l(ii,3),1)*phi_l+Nodes_rz(l(ii,4),1)*phi_c+Nodes_rz(l(ii,1),1)*phi_r;
    yvec = Nodes_rz(l(ii,3),2)*phi_l+Nodes_rz(l(ii,4),2)*phi_c+Nodes_rz(l(ii,1),2)*phi_r;
    plot(xvec,yvec,'-k')
end
%plotting element numbers
for ii = 1:n_el
    text(mean(Nodes_rz(l(ii,[4 5 6]),1))-.03,mean(Nodes_rz(l(ii,[4 5 6]),2)), ...
         ['e = ',num2str(ii)],'Fontsize',24)
end
%plotting the velocity-only nodes
for ii = 1:n_el
    plot(Nodes_rz(l(ii,[4 5 6]),1),Nodes_rz(l(ii,[4 5 6]),2),'sk', ...
         'MarkerSize',6)
end
%plotting the pressure-velocity nodes
for ii = 1:n_el
    plot(Nodes_rz(l(ii,[1 2 3]),1),Nodes_rz(l(ii,[1 2 3]),2),'ok', ...
        'MarkerSize',6)
end
%plotting global node numbers
for ii = 1:n_v
    text(Nodes_rz(ii,1),Nodes_rz(ii,2),[' i = ',num2str(ii)],'FontSize',16)
end
set(gca,'FontSize',24)
title('Mesh with the element numbers')












%element boundaries and node numbers
h3 = figure;
hold on
for ii = 1:n_el
    xvec = Nodes_rz(l(ii,1),1)*phi_l+Nodes_rz(l(ii,5),1)*phi_c+Nodes_rz(l(ii,2),1)*phi_r;
    yvec = Nodes_rz(l(ii,1),2)*phi_l+Nodes_rz(l(ii,5),2)*phi_c+Nodes_rz(l(ii,2),2)*phi_r;
    plot(xvec,yvec,'-k')
    xvec = Nodes_rz(l(ii,2),1)*phi_l+Nodes_rz(l(ii,6),1)*phi_c+Nodes_rz(l(ii,3),1)*phi_r;
    yvec = Nodes_rz(l(ii,2),2)*phi_l+Nodes_rz(l(ii,6),2)*phi_c+Nodes_rz(l(ii,3),2)*phi_r;
    plot(xvec,yvec,'-k')
    xvec = Nodes_rz(l(ii,3),1)*phi_l+Nodes_rz(l(ii,4),1)*phi_c+Nodes_rz(l(ii,1),1)*phi_r;
    yvec = Nodes_rz(l(ii,3),2)*phi_l+Nodes_rz(l(ii,4),2)*phi_c+Nodes_rz(l(ii,1),2)*phi_r;
    plot(xvec,yvec,'-k')
end
%plotting the velocity only nodes
for ii = 1:n_el
    plot(Nodes_rz(l(ii,[4 5 6]),1),Nodes_rz(l(ii,[4 5 6]),2),'sk','MarkerSize',6)
end
%plotting the pressure-velocity nodes
for ii = 1:n_el
    plot(Nodes_rz(l(ii,[1 2 3]),1),Nodes_rz(l(ii,[1 2 3]),2),'ok','MarkerSize',6)
end
%plotting global node numbers
for ii = 1:n_v
    text(Nodes_rz(ii,1),Nodes_rz(ii,2),[' i = ',num2str(ii)])
end
set(gca,'FontSize',24)
title('Mesh with node numbers')











%elements and local to globa pressure nodes map
h4 = figure;
hold on
for ii = 1:n_el
    xvec = Nodes_rz(l(ii,1),1)*phi_l+Nodes_rz(l(ii,5),1)*phi_c+Nodes_rz(l(ii,2),1)*phi_r;
    yvec = Nodes_rz(l(ii,1),2)*phi_l+Nodes_rz(l(ii,5),2)*phi_c+Nodes_rz(l(ii,2),2)*phi_r;
    plot(xvec,yvec,'-k')
    xvec = Nodes_rz(l(ii,2),1)*phi_l+Nodes_rz(l(ii,6),1)*phi_c+Nodes_rz(l(ii,3),1)*phi_r;
    yvec = Nodes_rz(l(ii,2),2)*phi_l+Nodes_rz(l(ii,6),2)*phi_c+Nodes_rz(l(ii,3),2)*phi_r;
    plot(xvec,yvec,'-k')
    xvec = Nodes_rz(l(ii,3),1)*phi_l+Nodes_rz(l(ii,4),1)*phi_c+Nodes_rz(l(ii,1),1)*phi_r;
    yvec = Nodes_rz(l(ii,3),2)*phi_l+Nodes_rz(l(ii,4),2)*phi_c+Nodes_rz(l(ii,1),2)*phi_r;
    plot(xvec,yvec,'-k')
end
%plotting the pressure-velocity nodes
tweak = .35;%parameter between 0 en 1 to choose location of text labels for 
%nodes in more than one element
for ii = 1:n_el
    plot(Nodes_rz(l(ii,:),1),Nodes_rz(l(ii,:),2),'ok','MarkerSize',6)
end
for ii = 1:n_el
    text((Nodes_rz(l(ii,1),1)+tweak*mean(Nodes_rz(l(ii,:),1)))/(1+tweak), ...
         (Nodes_rz(l(ii,1),2)+tweak*mean(Nodes_rz(l(ii,:),2)))/(1+tweak), ...
         ['l(',num2str(ii),',',num2str(1),') = ',num2str(l(ii,1))],'Fontsize',10);
    text((Nodes_rz(l(ii,2),1)+tweak*mean(Nodes_rz(l(ii,:),1)))/(1+tweak), ...
         (Nodes_rz(l(ii,2),2)+tweak*mean(Nodes_rz(l(ii,:),2)))/(1+tweak), ...
         ['l(',num2str(ii),',',num2str(2),') = ',num2str(l(ii,2))],'Fontsize',10);
    text((Nodes_rz(l(ii,3),1)+tweak*mean(Nodes_rz(l(ii,:),1)))/(1+tweak), ...
         (Nodes_rz(l(ii,3),2)+tweak*mean(Nodes_rz(l(ii,:),2)))/(1+tweak), ...
         ['l(',num2str(ii),',',num2str(3),') = ',num2str(l(ii,3))],'Fontsize',10);
    text((Nodes_rz(l(ii,4),1)+tweak*mean(Nodes_rz(l(ii,:),1)))/(1+tweak), ...
         (Nodes_rz(l(ii,4),2)+tweak*mean(Nodes_rz(l(ii,:),2)))/(1+tweak), ...
         ['l(',num2str(ii),',',num2str(4),') = ',num2str(l(ii,4))],'Fontsize',10);
    text((Nodes_rz(l(ii,5),1)+tweak*mean(Nodes_rz(l(ii,:),1)))/(1+tweak), ...
         (Nodes_rz(l(ii,5),2)+tweak*mean(Nodes_rz(l(ii,:),2)))/(1+tweak), ...
         ['l(',num2str(ii),',',num2str(5),') = ',num2str(l(ii,5))],'Fontsize',10);
    text((Nodes_rz(l(ii,6),1)+tweak*mean(Nodes_rz(l(ii,:),1)))/(1+tweak), ...
         (Nodes_rz(l(ii,6),2)+tweak*mean(Nodes_rz(l(ii,:),2)))/(1+tweak), ...
         ['l(',num2str(ii),',',num2str(6),') = ',num2str(l(ii,6))],'Fontsize',10);
end
%plotting element numbers
for ii = 1:n_el
    text(mean(Nodes_rz(l(ii,[4 5 6]),1))-.03,mean(Nodes_rz(l(ii,[4 5 6]),2)), ...
         ['e = ',num2str(ii)],'Fontsize',24)
end
set(gca,'FontSize',24)
title('Mesh with local to global velocity node number maps')












%elements and pressure-node numbers
h5 = figure;
hold on
for ii = 1:n_el
    xvec = Nodes_rz(l(ii,1),1)*phi_l+Nodes_rz(l(ii,5),1)*phi_c+Nodes_rz(l(ii,2),1)*phi_r;
    yvec = Nodes_rz(l(ii,1),2)*phi_l+Nodes_rz(l(ii,5),2)*phi_c+Nodes_rz(l(ii,2),2)*phi_r;
    plot(xvec,yvec,'-k')
    xvec = Nodes_rz(l(ii,2),1)*phi_l+Nodes_rz(l(ii,6),1)*phi_c+Nodes_rz(l(ii,3),1)*phi_r;
    yvec = Nodes_rz(l(ii,2),2)*phi_l+Nodes_rz(l(ii,6),2)*phi_c+Nodes_rz(l(ii,3),2)*phi_r;
    plot(xvec,yvec,'-k')
    xvec = Nodes_rz(l(ii,3),1)*phi_l+Nodes_rz(l(ii,4),1)*phi_c+Nodes_rz(l(ii,1),1)*phi_r;
    yvec = Nodes_rz(l(ii,3),2)*phi_l+Nodes_rz(l(ii,4),2)*phi_c+Nodes_rz(l(ii,1),2)*phi_r;
    plot(xvec,yvec,'-k')
end
%plotting the pressure-velocity nodes
tweakx = .01;
tweaky = .01;
for ii = 1:n_el
    plot(Nodes_rz(l(ii,1:3),1),Nodes_rz(l(ii,1:3),2),'ok','MarkerSize',6)
end
for ii = 1:n_el
    text(Nodes_rz(l(ii,1),1)+tweakx, Nodes_rz(l(ii,1),2)+tweaky, ['k = ',num2str(lp(ii,1))],'Fontsize',16);
    text(Nodes_rz(l(ii,2),1)+tweakx, Nodes_rz(l(ii,2),2)+tweaky, ['k = ',num2str(lp(ii,2))],'Fontsize',16);
    text(Nodes_rz(l(ii,3),1)+tweakx, Nodes_rz(l(ii,3),2)+tweaky, ['k = ',num2str(lp(ii,3))],'Fontsize',16);
end
%plotting element numbers
for ii = 1:n_el
    text(mean(Nodes_rz(l(ii,[4 5 6]),1))-.03,mean(Nodes_rz(l(ii,[4 5 6]),2)), ...
         ['e = ',num2str(ii)],'Fontsize',24)
end
set(gca,'Fontsize',24)
title('Mesh with local pressure-node number','Fontsize',24)











%elements and local to global pressure nodes map
h6 = figure;
hold on
for ii = 1:n_el
    xvec = Nodes_rz(l(ii,1),1)*phi_l+Nodes_rz(l(ii,5),1)*phi_c+Nodes_rz(l(ii,2),1)*phi_r;
    yvec = Nodes_rz(l(ii,1),2)*phi_l+Nodes_rz(l(ii,5),2)*phi_c+Nodes_rz(l(ii,2),2)*phi_r;
    plot(xvec,yvec,'-k')
    xvec = Nodes_rz(l(ii,2),1)*phi_l+Nodes_rz(l(ii,6),1)*phi_c+Nodes_rz(l(ii,3),1)*phi_r;
    yvec = Nodes_rz(l(ii,2),2)*phi_l+Nodes_rz(l(ii,6),2)*phi_c+Nodes_rz(l(ii,3),2)*phi_r;
    plot(xvec,yvec,'-k')
    xvec = Nodes_rz(l(ii,3),1)*phi_l+Nodes_rz(l(ii,4),1)*phi_c+Nodes_rz(l(ii,1),1)*phi_r;
    yvec = Nodes_rz(l(ii,3),2)*phi_l+Nodes_rz(l(ii,4),2)*phi_c+Nodes_rz(l(ii,1),2)*phi_r;
    plot(xvec,yvec,'-k')
end
%plotting the pressure-velocity nodes
tweak = .7;%parameter between 0 en 1 to choose location of text labels for 
%nodes in more than one element
for ii = 1:n_el
    plot(Nodes_rz(l(ii,1:3),1),Nodes_rz(l(ii,1:3),2),'ok','MarkerSize',6)
end
for ii = 1:n_el
    text((Nodes_rz(l(ii,1),1)+tweak*mean(Nodes_rz(l(ii,:),1)))/(1+tweak), ...
         (Nodes_rz(l(ii,1),2)+tweak*mean(Nodes_rz(l(ii,:),2)))/(1+tweak), ...
         ['lp(',num2str(ii),',',num2str(1),') = ',num2str(lp(ii,1))]);
    text((Nodes_rz(l(ii,2),1)+tweak*mean(Nodes_rz(l(ii,:),1)))/(1+tweak), ...
         (Nodes_rz(l(ii,2),2)+tweak*mean(Nodes_rz(l(ii,:),2)))/(1+tweak), ...
         ['lp(',num2str(ii),',',num2str(2),') = ',num2str(lp(ii,2))]);
    text((Nodes_rz(l(ii,3),1)+tweak*mean(Nodes_rz(l(ii,:),1)))/(1+tweak), ...
         (Nodes_rz(l(ii,3),2)+tweak*mean(Nodes_rz(l(ii,:),2)))/(1+tweak), ...
         ['lp(',num2str(ii),',',num2str(3),') = ',num2str(lp(ii,3))]);
end
% %plotting element numbers
for ii = 1:n_el
    text(mean(Nodes_rz(l(ii,[4 5 6]),1))-.03,mean(Nodes_rz(l(ii,[4 5 6]),2)), ...
         ['e = ',num2str(ii)],'Fontsize',24)
end
set(gca,'FontSize',24)
title('Mesh with local to global pressure node number maps')











%elements and local to global boundary 1 nodes map
h7 = figure;
hold on
for ii = 1:n_el
    xvec = Nodes_rz(l(ii,1),1)*phi_l+Nodes_rz(l(ii,5),1)*phi_c+Nodes_rz(l(ii,2),1)*phi_r;
    yvec = Nodes_rz(l(ii,1),2)*phi_l+Nodes_rz(l(ii,5),2)*phi_c+Nodes_rz(l(ii,2),2)*phi_r;
    plot(xvec,yvec,'-k')
    xvec = Nodes_rz(l(ii,2),1)*phi_l+Nodes_rz(l(ii,6),1)*phi_c+Nodes_rz(l(ii,3),1)*phi_r;
    yvec = Nodes_rz(l(ii,2),2)*phi_l+Nodes_rz(l(ii,6),2)*phi_c+Nodes_rz(l(ii,3),2)*phi_r;
    plot(xvec,yvec,'-k')
    xvec = Nodes_rz(l(ii,3),1)*phi_l+Nodes_rz(l(ii,4),1)*phi_c+Nodes_rz(l(ii,1),1)*phi_r;
    yvec = Nodes_rz(l(ii,3),2)*phi_l+Nodes_rz(l(ii,4),2)*phi_c+Nodes_rz(l(ii,1),2)*phi_r;
    plot(xvec,yvec,'-k')
end
%plotting the pressure-velocity nodes
tweak = .3;%parameter between 0 en 1 to choose location of text labels for 
%nodes in more than one element
shift = .05;
for ii = 1:n1_el
    plot(Nodes_rz(l_1(ii,:),1),Nodes_rz(l_1(ii,:),2),'ok','MarkerSize',6)
end
for ii = 1:n1_el
    text((Nodes_rz(l_1(ii,1),1)+tweak*mean(Nodes_rz(l_1(ii,:),1)))/(1+tweak), ...
         (Nodes_rz(l_1(ii,1),2)+tweak*mean(Nodes_rz(l_1(ii,:),2)))/(1+tweak)+shift, ...
         ['l_1(',num2str(ii),',',num2str(1),') = ',num2str(l_1(ii,1))],'rotation',90);
    text((Nodes_rz(l_1(ii,2),1)+tweak*mean(Nodes_rz(l_1(ii,:),1)))/(1+tweak), ...
         (Nodes_rz(l_1(ii,2),2)+tweak*mean(Nodes_rz(l_1(ii,:),2)))/(1+tweak)+shift, ...
         ['l_1(',num2str(ii),',',num2str(2),') = ',num2str(l_1(ii,2))],'rotation',90);
    text((Nodes_rz(l_1(ii,3),1)+tweak*mean(Nodes_rz(l_1(ii,:),1)))/(1+tweak), ...
         (Nodes_rz(l_1(ii,3),2)+tweak*mean(Nodes_rz(l_1(ii,:),2)))/(1+tweak)+shift, ...
         ['l_1(',num2str(ii),',',num2str(3),') = ',num2str(l_1(ii,3))],'rotation',90);
end
set(gca,'FontSize',24)
title('Mesh with local to global boundary 1 node number maps')






%elements and local to boundary 1 nodes map
h8 = figure;
hold on
for ii = 1:n_el
    xvec = Nodes_rz(l(ii,1),1)*phi_l+Nodes_rz(l(ii,5),1)*phi_c+Nodes_rz(l(ii,2),1)*phi_r;
    yvec = Nodes_rz(l(ii,1),2)*phi_l+Nodes_rz(l(ii,5),2)*phi_c+Nodes_rz(l(ii,2),2)*phi_r;
    plot(xvec,yvec,'-k')
    xvec = Nodes_rz(l(ii,2),1)*phi_l+Nodes_rz(l(ii,6),1)*phi_c+Nodes_rz(l(ii,3),1)*phi_r;
    yvec = Nodes_rz(l(ii,2),2)*phi_l+Nodes_rz(l(ii,6),2)*phi_c+Nodes_rz(l(ii,3),2)*phi_r;
    plot(xvec,yvec,'-k')
    xvec = Nodes_rz(l(ii,3),1)*phi_l+Nodes_rz(l(ii,4),1)*phi_c+Nodes_rz(l(ii,1),1)*phi_r;
    yvec = Nodes_rz(l(ii,3),2)*phi_l+Nodes_rz(l(ii,4),2)*phi_c+Nodes_rz(l(ii,1),2)*phi_r;
    plot(xvec,yvec,'-k')
end
%plotting the pressure-velocity nodes
tweak = .3;%parameter between 0 en 1 to choose location of text labels for 
%nodes in more than one element
shift = .05;
for ii = 1:n1_el
    plot(Nodes_rz(l_1(ii,:),1),Nodes_rz(l_1(ii,:),2),'ok','MarkerSize',6)
end
% axis equal
for ii = 1:n1_el
    text((Nodes_rz(l_1(ii,1),1)+tweak*mean(Nodes_rz(l_1(ii,:),1)))/(1+tweak), ...
         (Nodes_rz(l_1(ii,1),2)+tweak*mean(Nodes_rz(l_1(ii,:),2)))/(1+tweak)+shift, ...
         ['l^1_1(',num2str(ii),',',num2str(1),') = ',num2str(l1_1(ii,1))],'rotation',90);
    text((Nodes_rz(l_1(ii,2),1)+tweak*mean(Nodes_rz(l_1(ii,:),1)))/(1+tweak), ...
         (Nodes_rz(l_1(ii,2),2)+tweak*mean(Nodes_rz(l_1(ii,:),2)))/(1+tweak)+shift, ...
         ['l^1_1(',num2str(ii),',',num2str(2),') = ',num2str(l1_1(ii,2))],'rotation',90);
    text((Nodes_rz(l_1(ii,3),1)+tweak*mean(Nodes_rz(l_1(ii,:),1)))/(1+tweak), ...
         (Nodes_rz(l_1(ii,3),2)+tweak*mean(Nodes_rz(l_1(ii,:),2)))/(1+tweak)+shift, ...
         ['l^1_1(',num2str(ii),',',num2str(3),') = ',num2str(l1_1(ii,3))],'rotation',90);
end
set(gca,'FontSize',24)%,'ylim',[0 1],'xlim',[-1 0.2])
% grid on
% axis off
title('Mesh with local to boundary-1-node number map')





%elements and local to global nodes map
h9 = figure;
hold on
for ii = 1:n_el
    xvec = Nodes_rz(l(ii,1),1)*phi_l+Nodes_rz(l(ii,5),1)*phi_c+Nodes_rz(l(ii,2),1)*phi_r;
    yvec = Nodes_rz(l(ii,1),2)*phi_l+Nodes_rz(l(ii,5),2)*phi_c+Nodes_rz(l(ii,2),2)*phi_r;
    plot(xvec,yvec,'-k')
    xvec = Nodes_rz(l(ii,2),1)*phi_l+Nodes_rz(l(ii,6),1)*phi_c+Nodes_rz(l(ii,3),1)*phi_r;
    yvec = Nodes_rz(l(ii,2),2)*phi_l+Nodes_rz(l(ii,6),2)*phi_c+Nodes_rz(l(ii,3),2)*phi_r;
    plot(xvec,yvec,'-k')
    xvec = Nodes_rz(l(ii,3),1)*phi_l+Nodes_rz(l(ii,4),1)*phi_c+Nodes_rz(l(ii,1),1)*phi_r;
    yvec = Nodes_rz(l(ii,3),2)*phi_l+Nodes_rz(l(ii,4),2)*phi_c+Nodes_rz(l(ii,1),2)*phi_r;
    plot(xvec,yvec,'-k')
end
%plotting the pressure-velocity nodes
tweak = .3;%parameter between 0 en 1 to choose location of text labels for 
% %nodes in more than one element
shift = .05;
for ii = 1:n2_el
    plot(Nodes_rz(l_2(ii,:),1),Nodes_rz(l_2(ii,:),2),'ok','MarkerSize',6)
end
for ii = 1:n2_el
    text((Nodes_rz(l_2(ii,1),1)+tweak*mean(Nodes_rz(l_2(ii,:),1)))/(1+tweak), ...
          Nodes_rz(l_2(ii,1),2)+shift, ...
         ['l_2(',num2str(ii),',',num2str(1),') = ',num2str(l_2(ii,1))],'rotation',90);
    text((Nodes_rz(l_2(ii,2),1)+tweak*mean(Nodes_rz(l_2(ii,:),1)))/(1+tweak), ...
          Nodes_rz(l_2(ii,2),2)+shift, ...
         ['l_2(',num2str(ii),',',num2str(2),') = ',num2str(l_2(ii,2))],'rotation',90);
    text((Nodes_rz(l_2(ii,3),1)+tweak*mean(Nodes_rz(l_2(ii,:),1)))/(1+tweak), ...
          Nodes_rz(l_2(ii,3),2)+shift, ...
         ['l_2(',num2str(ii),',',num2str(3),') = ',num2str(l_2(ii,3))],'rotation',90);
end
set(gca,'FontSize',24)%,'ylim',[0 1],'xlim',[-1 0.2])
title('Mesh with local to global boundary-2-node number maps')





%elements and local to global boundary 2 nodes map
h10 = figure;
hold on
for ii = 1:n_el
    xvec = Nodes_rz(l(ii,1),1)*phi_l+Nodes_rz(l(ii,5),1)*phi_c+Nodes_rz(l(ii,2),1)*phi_r;
    yvec = Nodes_rz(l(ii,1),2)*phi_l+Nodes_rz(l(ii,5),2)*phi_c+Nodes_rz(l(ii,2),2)*phi_r;
    plot(xvec,yvec,'-k')
    xvec = Nodes_rz(l(ii,2),1)*phi_l+Nodes_rz(l(ii,6),1)*phi_c+Nodes_rz(l(ii,3),1)*phi_r;
    yvec = Nodes_rz(l(ii,2),2)*phi_l+Nodes_rz(l(ii,6),2)*phi_c+Nodes_rz(l(ii,3),2)*phi_r;
    plot(xvec,yvec,'-k')
    xvec = Nodes_rz(l(ii,3),1)*phi_l+Nodes_rz(l(ii,4),1)*phi_c+Nodes_rz(l(ii,1),1)*phi_r;
    yvec = Nodes_rz(l(ii,3),2)*phi_l+Nodes_rz(l(ii,4),2)*phi_c+Nodes_rz(l(ii,1),2)*phi_r;
    plot(xvec,yvec,'-k')
end
%plotting the pressure-velocity nodes
tweak = .3;%parameter between 0 en 1 to choose location of text labels for 
% %nodes in more than one element
shift = .05;
for ii = 1:n2_el
    plot(Nodes_rz(l_2(ii,:),1),Nodes_rz(l_2(ii,:),2),'ok','MarkerSize',6)
end
for ii = 1:n2_el
    text((Nodes_rz(l_2(ii,1),1)+tweak*mean(Nodes_rz(l_2(ii,:),1)))/(1+tweak), ...
          Nodes_rz(l_2(ii,1),2)+shift, ...
         ['l^2_2(',num2str(ii),',',num2str(1),') = ',num2str(l2_2(ii,1))],'rotation',90);
    text((Nodes_rz(l_2(ii,2),1)+tweak*mean(Nodes_rz(l_2(ii,:),1)))/(1+tweak), ...
          Nodes_rz(l_2(ii,2),2)+shift, ...
         ['l^2_2(',num2str(ii),',',num2str(2),') = ',num2str(l2_2(ii,2))],'rotation',90);
    text((Nodes_rz(l_2(ii,3),1)+tweak*mean(Nodes_rz(l_2(ii,:),1)))/(1+tweak), ...
          Nodes_rz(l_2(ii,3),2)+shift, ...
         ['l^2_2(',num2str(ii),',',num2str(3),') = ',num2str(l2_2(ii,3))],'rotation',90);
end
set(gca,'FontSize',24)
title('Mesh with local to boundary-2-node number maps')







%boundary 3 elements and local to global nodes map
h11 = figure;
hold on
for ii = 1:n_el
    xvec = Nodes_rz(l(ii,1),1)*phi_l+Nodes_rz(l(ii,5),1)*phi_c+Nodes_rz(l(ii,2),1)*phi_r;
    yvec = Nodes_rz(l(ii,1),2)*phi_l+Nodes_rz(l(ii,5),2)*phi_c+Nodes_rz(l(ii,2),2)*phi_r;
    plot(xvec,yvec,'-k')
    xvec = Nodes_rz(l(ii,2),1)*phi_l+Nodes_rz(l(ii,6),1)*phi_c+Nodes_rz(l(ii,3),1)*phi_r;
    yvec = Nodes_rz(l(ii,2),2)*phi_l+Nodes_rz(l(ii,6),2)*phi_c+Nodes_rz(l(ii,3),2)*phi_r;
    plot(xvec,yvec,'-k')
    xvec = Nodes_rz(l(ii,3),1)*phi_l+Nodes_rz(l(ii,4),1)*phi_c+Nodes_rz(l(ii,1),1)*phi_r;
    yvec = Nodes_rz(l(ii,3),2)*phi_l+Nodes_rz(l(ii,4),2)*phi_c+Nodes_rz(l(ii,1),2)*phi_r;
    plot(xvec,yvec,'-k')
end
%plotting the pressure-velocity nodes
tweak = .3;%parameter between 0 en 1 to choose location of text labels for 
% %nodes in more than one element
shift = .05;
for ii = 1:n3_el
    plot(Nodes_rz(l_3(ii,:),1),Nodes_rz(l_3(ii,:),2),'ok','MarkerSize',6)
end
for ii = 1:n3_el
    text(Nodes_rz(l_3(ii,1),1)+shift, ...
         (Nodes_rz(l_3(ii,1),2)+tweak*mean(Nodes_rz(l_3(ii,:),2)))/(1+tweak), ...
         ['l_3(',num2str(ii),',',num2str(1),') = ',num2str(l_3(ii,1))]);
    text(Nodes_rz(l_3(ii,2),1)+shift, ...
         (Nodes_rz(l_3(ii,2),2)+tweak*mean(Nodes_rz(l_3(ii,:),2)))/(1+tweak), ...
         ['l_3(',num2str(ii),',',num2str(2),') = ',num2str(l_3(ii,2))]);
    text(Nodes_rz(l_3(ii,3),1)+shift, ...
         (Nodes_rz(l_3(ii,3),2)+tweak*mean(Nodes_rz(l_3(ii,:),2)))/(1+tweak), ...
         ['l_3(',num2str(ii),',',num2str(3),') = ',num2str(l_3(ii,3))]);
end
set(gca,'FontSize',24)%,'ylim',[0 1],'xlim',[-1 0.2])
title('Mesh with local to global boundary-3-node number maps')





%elements and local to global boundary 3 nodes map
h12 = figure;
hold on
for ii = 1:n_el
    xvec = Nodes_rz(l(ii,1),1)*phi_l+Nodes_rz(l(ii,5),1)*phi_c+Nodes_rz(l(ii,2),1)*phi_r;
    yvec = Nodes_rz(l(ii,1),2)*phi_l+Nodes_rz(l(ii,5),2)*phi_c+Nodes_rz(l(ii,2),2)*phi_r;
    plot(xvec,yvec,'-k')
    xvec = Nodes_rz(l(ii,2),1)*phi_l+Nodes_rz(l(ii,6),1)*phi_c+Nodes_rz(l(ii,3),1)*phi_r;
    yvec = Nodes_rz(l(ii,2),2)*phi_l+Nodes_rz(l(ii,6),2)*phi_c+Nodes_rz(l(ii,3),2)*phi_r;
    plot(xvec,yvec,'-k')
    xvec = Nodes_rz(l(ii,3),1)*phi_l+Nodes_rz(l(ii,4),1)*phi_c+Nodes_rz(l(ii,1),1)*phi_r;
    yvec = Nodes_rz(l(ii,3),2)*phi_l+Nodes_rz(l(ii,4),2)*phi_c+Nodes_rz(l(ii,1),2)*phi_r;
    plot(xvec,yvec,'-k')
end
%plotting the pressure-velocity nodes
tweak = .3;%parameter between 0 en 1 to choose location of text labels for 
% %nodes in more than one element
shift = .05;
for ii = 1:n3_el
    plot(Nodes_rz(l_3(ii,:),1),Nodes_rz(l_3(ii,:),2),'ok','MarkerSize',6)
end
for ii = 1:n3_el
    text(Nodes_rz(l_3(ii,1),1)+shift, ...
         (Nodes_rz(l_3(ii,1),2)+tweak*mean(Nodes_rz(l_3(ii,:),2)))/(1+tweak), ...
         ['l^3_3(',num2str(ii),',',num2str(1),') = ',num2str(l3_3(ii,1))]);
    text(Nodes_rz(l_3(ii,2),1)+shift, ...
         (Nodes_rz(l_3(ii,2),2)+tweak*mean(Nodes_rz(l_3(ii,:),2)))/(1+tweak), ...
         ['l^3_3(',num2str(ii),',',num2str(2),') = ',num2str(l3_3(ii,2))]);
    text(Nodes_rz(l_3(ii,3),1)+shift, ...
         (Nodes_rz(l_3(ii,3),2)+tweak*mean(Nodes_rz(l_3(ii,:),2)))/(1+tweak), ...
         ['l^3_3(',num2str(ii),',',num2str(3),') = ',num2str(l3_3(ii,3))]);
end
set(gca,'FontSize',24)
title('Mesh with local to boundary-3-node number maps')





%effect of perturbing f
h13 = figure;
hold on
delt_s = .02;
%plotting mesh of elements
for ii = 1:n_el
    xvec = Nodes_rz(l(ii,1),1)*phi_l+Nodes_rz(l(ii,5),1)*phi_c+Nodes_rz(l(ii,2),1)*phi_r;
    yvec = Nodes_rz(l(ii,1),2)*phi_l+Nodes_rz(l(ii,5),2)*phi_c+Nodes_rz(l(ii,2),2)*phi_r;
    plot(xvec,yvec,'-','color',[.7 .7 .7],'Linewidth',4)
    xvec = Nodes_rz(l(ii,2),1)*phi_l+Nodes_rz(l(ii,6),1)*phi_c+Nodes_rz(l(ii,3),1)*phi_r;
    yvec = Nodes_rz(l(ii,2),2)*phi_l+Nodes_rz(l(ii,6),2)*phi_c+Nodes_rz(l(ii,3),2)*phi_r;
    plot(xvec,yvec,'-','color',[.7 .7 .7],'Linewidth',4)
    xvec = Nodes_rz(l(ii,3),1)*phi_l+Nodes_rz(l(ii,4),1)*phi_c+Nodes_rz(l(ii,1),1)*phi_r;
    yvec = Nodes_rz(l(ii,3),2)*phi_l+Nodes_rz(l(ii,4),2)*phi_c+Nodes_rz(l(ii,1),2)*phi_r;
    plot(xvec,yvec,'-','color',[.7 .7 .7],'Linewidth',4)
end
%plotting element numbers
for ii = 1:n_el
    text(mean(Nodes_rz(l(ii,[4 5 6]),1))-.03,mean(Nodes_rz(l(ii,[4 5 6]),2)), ...
         ['e = ',num2str(ii)],'Fontsize',12)
end
%plotting the velocity-only nodes
for ii = 1:n_el
    plot(Nodes_rz(l(ii,[4 5 6]),1),Nodes_rz(l(ii,[4 5 6]),2),'sk', ...
         'MarkerSize',6)
end
%plotting the pressure-velocity nodes
for ii = 1:n_el
    plot(Nodes_rz(l(ii,[1 2 3]),1),Nodes_rz(l(ii,[1 2 3]),2),'ok', ...
        'MarkerSize',6)
end
%plotting global node numbers
for ii = 1:n_v
    text(Nodes_rz(ii,1),Nodes_rz(ii,2),[' i = ',num2str(ii)],'FontSize',12)
end
new_nodes = ...
    nodes_relocator_split_v02([spine_lengths_ini(1)*(1+delt_s); ...
                                       spine_lengths_ini(2:end)], ...
                                       spine_feet_locator,n_spines, ...
                                       n_spines_incr,n_v,alpha_s, ...
                                       n_nodes_const_per_spine, ...
                                       n_el_sing);
for ii = 1:n_el
    xvec = new_nodes(l(ii,1),1)*phi_l+new_nodes(l(ii,5),1)*phi_c+new_nodes(l(ii,2),1)*phi_r;
    yvec = new_nodes(l(ii,1),2)*phi_l+new_nodes(l(ii,5),2)*phi_c+new_nodes(l(ii,2),2)*phi_r;
    plot(xvec,yvec,'-','color',[0 0 0],'Linewidth',1)
    xvec = new_nodes(l(ii,2),1)*phi_l+new_nodes(l(ii,6),1)*phi_c+new_nodes(l(ii,3),1)*phi_r;
    yvec = new_nodes(l(ii,2),2)*phi_l+new_nodes(l(ii,6),2)*phi_c+new_nodes(l(ii,3),2)*phi_r;
    plot(xvec,yvec,'-','color',[0 0 0],'Linewidth',1)
    xvec = new_nodes(l(ii,3),1)*phi_l+new_nodes(l(ii,4),1)*phi_c+new_nodes(l(ii,1),1)*phi_r;
    yvec = new_nodes(l(ii,3),2)*phi_l+new_nodes(l(ii,4),2)*phi_c+new_nodes(l(ii,1),2)*phi_r;
    plot(xvec,yvec,'-','color',[0 0 0],'Linewidth',1)
end
set(gca,'FontSize',24,'xlim',[-2.2 .2],'ylim',[0 1.05])
title('Effect of perturbing $f$','interpreter','Latex')



