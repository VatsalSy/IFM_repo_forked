close all
clear
clc

load('bound_to_glob2_map.mat')
load('it.mat')
load('spine_feet_locator.mat')
load('n_spines.mat')
load('n_spines_incr.mat')
load('n_v.mat')
load('Ts.mat')
load('rhos2_0.mat')
rhos2_ini = rhos2(1);
load('Ds.mat')

density = figure;
for ii = 0:it
    load(['spine_lengths_',num2str(ii),'.mat'])
    [Nodes_rz,~] ...
        = nodes_relocator_split_v01(spine_lengths, ...
                                    spine_feet_locator, ...
                                    n_spines, ...
                                    n_spines_incr,n_v);
    lim_left = Nodes_rz(1,1);
    lim_right = abs(Nodes_rz(bound_to_glob2_map(2),1)-Nodes_rz(1,1));
    semilogx([lim_left lim_right]/Ts,[rhos2_ini rhos2_ini])
    hold on
    semilogx([lim_left lim_right]/Ts,[Ds Ds])
    load(['rhos2_',num2str(ii),'.mat'])
    semilogx(flipud(abs(Nodes_rz(bound_to_glob2_map(2:end),1)-Nodes_rz(1,1)))/Ts, ...
             fliplr(rhos2(2:end)))
    grid on
    pause
    hold off
end
    