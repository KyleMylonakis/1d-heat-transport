clear; clc;
import Sample.*

%% Load in Processed data
n_systems = 11;
base_n = 16;
systems(1:n_systems)= Sample;
for j = 1:n_systems
   load(strcat(num2str(base_n + (j-1)*100),'_observables.mat')); 
   systems(j) = data;
   clear data;
end

%% Examine Processed data

for j = 1:n_systems
   figure
   subplot(3,2,1)
   p = plot(systems(j).m_loc_temp);
   title('Kinetic Temperatore')
   xlabel('Particle Number')
   ylabel('<p^2>')
   xlim([0 systems(j).m_n_particles])
   p.LineStyle = 'none';
   p.Marker = '.';
   p.MarkerSize = 4;
   
   subplot(3,2,2)
   p = plot(systems(j).m_loc_temp_no_loc_drift);
   title('Kinetic Temperatore - No local drift')
   xlabel('Particle Number')
   %ylabel('<p^2>')
   xlim([0 systems(j).m_n_particles])
   p.LineStyle = 'none';
   p.Marker = '.';
   p.MarkerSize = 4;
   
   subplot(3,2,3)
   p = plot(systems(j).m_loc_temp_no_tot_drift);
   title('Kinetic Temperatore -No total drift')
   xlabel('Particle Number')
   %ylabel('<p^2>')
   xlim([0 systems(j).m_n_particles])
   p.LineStyle = 'none';
   p.Marker = '.';
   p.MarkerSize = 4;
   
   subplot(3,2,4)
   fit = polyfit(1:(systems(j).m_n_particles-1),systems(j).m_loc_heat,1);
   hold on;
   p = plot(systems(j).m_loc_heat);
   q = plot(polyval(fit,1:(systems(j).m_n_particles-1)));
   hold off;
   title('Local Heat Flux')
   xlabel('Particle Number')
   %ylabel('<p^2>')
   xlim([0 systems(j).m_n_particles])
   p.LineStyle = 'none';
   p.Marker = '.';
   p.MarkerFaceColor = 'blue';
   p.MarkerSize = 6;
   q.LineWidth = 2;
   tx = text(systems(j).m_n_particles/2,max(systems(j).m_loc_heat),['Slope of linear best fit ' num2str(fit(1))]);
   tx.HorizontalAlignment = 'center';    % set horizontal alignment to center
   tx.VerticalAlignment = 'top';
   
   subplot(3,2,6)
   fit = polyfit(1:(systems(j).m_n_particles-1),systems(j).m_loc_heat_neglect_fluct,1);
   hold on;
   p = plot(systems(j).m_loc_heat_neglect_fluct);
   q = plot(polyval(fit,1:(systems(j).m_n_particles-1)));
   hold off;
   title('Local Heat Flux - Neglect local fluctuations')
   xlabel('Particle Number')
   %ylabel('<p^2>')
   xlim([0 systems(j).m_n_particles])
   p.LineStyle = 'none';
   p.Marker = '.';
   p.MarkerSize = 6;
   q.LineWidth = 2;
   
end