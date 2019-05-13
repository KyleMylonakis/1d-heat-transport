%% Define samples class

classdef Sample
   properties 
       a0 = 2^(1/6)
       m_n_particles
       m_disp
       m_vel
       m_pos
       m_num_samples
       m_inv_trunc_param
   end
   
   % Calculated Observables
   properties 
       m_loc_temp
       m_loc_temp_no_loc_drift
       m_loc_temp_no_tot_drift
       m_loc_energy
       m_loc_heat
       m_loc_heat_neglect_fluct
       m_tot_heat
       m_tot_heat_neglect_fluct
       m_conductivity
       m_conductivity_neglect_fluct
   end
   
   methods
       % Constructor
       function obj = Sample(disp, vel) 
           if nargin == 0
               fprintf('Using default constructor \n');
           else
               csv_size = size(vel);
               obj.m_n_particles = csv_size(2) - 1;
               obj.m_num_samples = csv_size(1);
               
               fprintf('Loading %i particle dispitions \n', obj.m_n_particles);
               obj.m_disp = disp(:,1:(end - 1));
               
               fprintf('Loading %i particle velocities \n', obj.m_n_particles);
               obj.m_vel = vel(:,1:(end - 1));
               
               fprintf('Constructing Positions from Displacements\n');
               for j = 1:obj.m_n_particles
                   obj.m_pos(:,j) = obj.m_disp(:,j) + obj.a0.*(j-1).*ones(obj.m_num_samples,1);
               end
               
               fprintf('Constructung local temp\n');
               obj.m_loc_temp = init_m_loc_temp(obj);
               
               fprintf('Constructung local temp with no local drift\n');
               obj.m_loc_temp_no_loc_drift = init_m_loc_temp_no_loc_drift(obj);
               
               fprintf('Constructing local temp with no total drift\n');
               obj.m_loc_temp_no_tot_drift = init_m_loc_temp_no_tot_drift(obj);
               
               fprintf('Constructing local energy\n');
               obj.m_loc_energy = init_m_loc_energy(obj);
               
               fprintf('Constructing local heat\n');
               obj.m_loc_heat = init_m_loc_heat(obj);
               
               fprintf('Constructing local heat ignoring local fluctuations\n');
               obj.m_loc_heat_neglect_fluct = init_m_loc_heat_neglect_fluct(obj);
               
               fprintf('Calculating total heat flux\n');
               obj.m_tot_heat = get_tot_heat(obj);
               
               fprintf('Calculating total heat flux ingnoring local fluctuations\n');
               obj.m_tot_heat_neglect_fluct = get_tot_heat_neglect_fluct(obj);
               
               fprintf('Calculating conductivity\n');
               obj.m_conductivity = get_conductivity(obj);
               
               fprintf('Calculating conductivity ignoring local fluctuations\n');
               obj.m_conductivity_neglect_fluct = get_conductivity_neglect_fluct(obj);
               
           end
       end

       % Calculated Observables
       function out = init_m_loc_temp(obj)
           out = mean(obj.m_vel.*obj.m_vel);
       end
       
       function out = get.m_loc_temp(obj)
           out = obj.m_loc_temp;
       end
       
       
       function out = init_m_loc_temp_no_loc_drift(obj)
           loc_drift = mean(obj.m_vel);
           out = mean( (obj.m_vel - repmat(loc_drift,obj.m_num_samples,1)).*(obj.m_vel - repmat(loc_drift,obj.m_num_samples,1) ) );
       end
       
       function out = get.m_loc_temp_no_loc_drift(obj)
           out = obj.m_loc_temp_no_loc_drift;
       end
       
       function out = init_m_loc_temp_no_tot_drift(obj)
           tot_drift = mean(mean(obj.m_vel));
           out = mean( (obj.m_vel - tot_drift.*ones(obj.m_num_samples,obj.m_n_particles)).*(obj.m_vel - tot_drift.*ones(obj.m_num_samples,obj.m_n_particles)));
       end
       
       function out = get.m_loc_temp_no_tot_drift(obj)
           out = obj.m_loc_temp_no_tot_drift;
       end
       
       function out = F(obj, x)
          out = -((-12 .*(x.^(-13))) + (6.*(x.^(-7))) + (4.*((x- obj.a0).^3)) ) ;
       end
       
       function out = V(obj, x)
           out = x.^(-12) - x.^(-6) + (x-obj.a0).^4;
       end
       
       function out = init_m_loc_energy(obj)
           k = (obj.m_vel.^2)/2;
           u = zeros(size(k));
           for j = 1:(obj.m_n_particles)
               if j == 1
                   u(:,1) = 0.5.*obj.V( obj.m_pos(:,2) - obj.m_pos(:,1) );
               elseif j < obj.m_n_particles
                   u(:,j) = 0.5.*(obj.V(obj.m_pos(:,j+1) - obj.m_pos(:,j)) ...
                         + obj.V( obj.m_pos(:,j) - obj.m_pos(:,j-1)));
               else
                   u(:,end) = 0.5.*obj.V(obj.m_pos(:,j) - obj.m_pos(:,j-1));
               end
           end
           out = k + u;
       end
       
       function out = get.m_loc_energy(obj)
           out = obj.m_loc_energy;
       end
       
       function out = init_m_loc_heat(obj)
           out = zeros(obj.m_num_samples,obj.m_n_particles-1);
           for j = 1:(obj.m_n_particles -1)
               out(:,j) = (0.5.*(obj.m_pos(:,j+1) - obj.m_pos(:,j)) ...
                   .*(obj.m_vel(:,j+1) + obj.m_vel(:,j)) ...
                   .*obj.F(obj.m_pos(:,j+1) - obj.m_pos(:,j)) ) ... 
                   + obj.m_vel(:,j).*obj.m_loc_energy(:,j);
           end
           out = mean(out);
       end
           
       function out = get.m_loc_heat(obj)
           out = obj.m_loc_heat;
       end
       
       function out = init_m_loc_heat_neglect_fluct(obj)
           out = zeros(obj.m_num_samples,obj.m_n_particles-1);
           for j = 1:(obj.m_n_particles -1)
               out(:,j) = (0.5.*(obj.a0) ...
                   .*(obj.m_vel(:,j+1) + obj.m_vel(:,j)) ...
                   .*obj.F(obj.m_pos(:,j+1) - obj.m_pos(:,j)) );
           end
           out = mean(out);
       end
       
       function obj = prepare_to_save(obj)
           fprintf('Removing non-averaged observables\n');
           obj.m_disp = NaN;
           obj.m_vel = NaN;
           obj.m_pos = NaN;
           obj.m_loc_energy = NaN;
       end
       
       function out = get_tot_heat(obj)
          out =  sum(obj.m_loc_heat);
       end
       
       function out = get_tot_heat_neglect_fluct(obj)
           out = sum(obj.m_loc_heat_neglect_fluct);
       end
       
       function out = get_conductivity(obj)
           out = -obj.m_tot_heat .* obj.m_n_particles ./ (obj.m_loc_temp(end) - obj.m_loc_temp(1));
       end
       
       function out = get_conductivity_neglect_fluct(obj)
           out = -obj.m_tot_heat_neglect_fluct .* obj.m_n_particles ./ (obj.m_loc_temp(end) - obj.m_loc_temp(1));
       end
   
       
   end
end