%% Reinitialize Variables
clear; clc;
import Sample.*

%% Import Data
% fprintf('Loading Data \n');
n_systems = 2;
base_n_part = 16;
for j = 2:n_systems %final j should be 11
    n_part = num2str(base_n_part + (j-1)*100);
    pos_file = strcat( './', n_part, '/pos_samp.csv');
    vel_file = strcat( './', n_part, '/vel_samp.csv');
    fprintf('Loading Data\n');
    data = Sample(base_n_part + (j-1)*100, ...
            csvread(pos_file), csvread(vel_file)).prepare_to_save();
    save(strcat(n_part, '_observables'), 'data');
    
end
    
%  samples = [
%  Sample(16, csvread('./16/pos_samp.csv'), csvread('./16/vel_samp.csv') )
% % Sample(116, csvread('./116/pos_samp.csv'), csvread('./116/vel_samp.csv') )
% % Sample(216, csvread('./216/pos_samp.csv'), csvread('./216/vel_samp.csv') )
% % Sample(316, csvread('./316/pos_samp.csv'), csvread('./316/vel_samp.csv') )
% % Sample(416, csvread('./416/pos_samp.csv'), csvread('./416/vel_samp.csv') )
% % Sample(516, csvread('./516/pos_samp.csv'), csvread('./516/vel_samp.csv') )
% % Sample(616, csvread('./616/pos_samp.csv'), csvread('./616/vel_samp.csv') )
% % Sample(716, csvread('./716/pos_samp.csv'), csvread('./716/vel_samp.csv') )
% % Sample(816, csvread('./816/pos_samp.csv'), csvread('./816/vel_samp.csv') )
% % Sample(916, csvread('./916/pos_samp.csv'), csvread('./916/vel_samp.csv') )
% % Sample(1016, csvread('./1016/pos_samp.csv'), csvread('./1016/vel_samp.csv')) 
% 
%  ];
% 
% fprintf('Data Loaded \n');

%% Compute Observables

% %a = samples(1).m_loc_temp
% %a = samples(1).m_loc_temp_no_loc_drift
% a = samples(1).m_loc_heat;
% b = samples(1).init_m_loc_heat_neglect_fluct;
% plot(abs(a-b)/(max(max(abs(a)))));
% samples(1) = samples(1).prepare_to_save();



