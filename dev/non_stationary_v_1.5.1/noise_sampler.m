%total_samples = 10000;
k_sqrd = 156*(2)^(-7.0/3.0) - 42*2^(-4.0/3.0);
total_samples = 1;
particles = 54000;
%particles = 500;
Ltemp = 0.2;
Rtemp = 0.8;
cutoff = particles - 100; % Ideally Number of particles needed to hit machine precision when evaluating 
mean = zeros(particles,1);

% Sample left bath positions
'Sampling Left Bath Positions'
left_pos_covariance = full(gallery('tridiag',particles,-1*k_sqrd,2*k_sqrd,-1*k_sqrd)); 
%left_pos_covariance = left_pos_covariance.*a0;
left_pos_covariance = Ltemp.*inv(left_pos_covariance);
left_pos_sample = mvnrnd(mean, left_pos_covariance, total_samples);
left_positions = left_pos_sample(:,1:cutoff)';
dlmwrite('left_noise_pos.dat', left_positions, 'precision', 16);
clear left_pos_covariance left_pos_sample left_positions;

% Sample right bath positions
'Sampling Right Bath Positions'
right_pos_covariance = full(gallery('tridiag',particles,-1*k_sqrd,2*k_sqrd,-1*k_sqrd)); 
%right_pos_covariance = right_pos_covariance.*a0;
right_pos_covariance = Rtemp.*inv(right_pos_covariance);
right_pos_sample = mvnrnd(mean,right_pos_covariance, total_samples);
right_positions = right_pos_sample(:,1:cutoff)';
dlmwrite('right_noise_pos.dat', right_positions, 'precision', 16);
clear right_pos_covariance right_pos_sample right_positions;

%left_pos_covariance = full(gallery('tridiag',particles,-1,2,-1)); 
%right_pos_covariance = full(gallery('tridiag',particles,-1,2,-1)); 
%left_pos_covariance = left_pos_covariance.*a0;
%right_pos_covariance = right_pos_covariance.*a0;

%'Inverting Covariance Matricies'
%left_pos_covariance = Ltemp.*inv(left_pos_covariance);
%right_pos_covariance = RTemp.*inv(right_pos_covariance);

%left_vel_covariance = 0.02.*ones(particles,particles); %kT = 0.01
%right_vel_covariance = 0.08.*ones(particles,particles); %kT = 0.01

%mean = zeros(particles,1);
%'Sampling Positions'
%left_pos_sample = mvnrnd(mean, left_pos_covariance, total_samples);
%right_pos_sample = mvnrnd(mean,right_pos_covariance, total_samples);
%left_vel_sample = mvnrnd(mean, left_vel_covariance, total_samples);
%temp=(mvnrnd(mean, left_vel_covariance, total_samples));
%var(temp)
'Sampling Velocities'
left_vel_sample =  normrnd(0, sqrt(Ltemp), [total_samples,particles]);
%var(left_vel_sample)
right_vel_sample = normrnd(0, sqrt(Rtemp), [total_samples,particles]); 
%right_vel_sample = mvnrnd(mean, right_vel_covariance, total_samples);



%left_positions = left_pos_sample(:,1:cutoff)';
%right_positions = right_pos_sample(:,1:cutoff)';
left_velocities = left_vel_sample(:,1:cutoff)';
right_velocities = right_vel_sample(:,1:cutoff)';

%save('sampled_noise.mat','left_positions', 'right_positions', 'left_velocities', 'right_velocities')
%dlmwrite('left_noise_pos.dat', left_positions, 'precision', 16);
dlmwrite('left_noise_vel.dat', left_velocities, 'precision', 16);
%dlmwrite('right_noise_pos.dat', right_positions, 'precision', 16);
dlmwrite('right_noise_vel.dat', right_velocities, 'precision', 16);
