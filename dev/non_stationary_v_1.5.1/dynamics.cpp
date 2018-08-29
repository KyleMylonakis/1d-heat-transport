/*
* Main file for the integrator of the system of ODE's
* First we import the random samples from the baths
* Then we use velocity Verlet and RK3 to simulate the dynamics
* After the time reaches a predefined equilibriation time, we sample the
* position and velocities and save the results to disk.
*/

#include "main.h"

int main()
{
    //std::cout << "dt_inv has value " << dt_inv << std::endl;


    if (((t_max -t_equi)*dt_inv) % sample_freq != 0)
    {
        std::cerr << "ERROR: LAST SAMPLES WILL CAUSE OVERFLOW" << std::endl;
        std::cerr << "CHECK THAT sample_freq DIVIDES (t_max -t_equi)*dt_inv " << std::endl;
        std::cerr << "EXITING" << std::endl;
        return 0;
    }

    if (dt_inv % 10 != 0)
    {
        std::cerr << "ERROR: INCORRECT SAMPLE FREQUENCE" << std::endl
            << "CHECK THAT 10 DIVIDES dt_inv" << std::endl;
        std::cerr << "EXITING" << std::endl;
        return 0;
    }

    // Load the randomly sampled initial positions and velocities for the
    // right and left baths
    std::cout << std::setprecision(std::numeric_limits<double>::max_digits10);
    std::cout << "Loading bath data" << std::endl;
    std::string line1, line2, line3, line4;
    std::ifstream left_noise_pos_f ("left_noise_pos.dat");
    std::ifstream left_noise_vel_f ("left_noise_vel.dat");
    std::ifstream right_noise_pos_f ("right_noise_pos.dat");
    std::ifstream right_noise_vel_f ("right_noise_vel.dat");

    if ( left_noise_pos_f.is_open() && left_noise_vel_f.is_open() 
            && right_noise_pos_f.is_open() && right_noise_vel_f.is_open() 
       )
    {
        for (long long int j = 0; j < n_bath_particles; ++j)
        {
            //std::cout << j << std::endl;
            std::getline(left_noise_pos_f,line1);
            std::getline(left_noise_vel_f,line2);
            std::getline(right_noise_pos_f,line3);
            std::getline(right_noise_vel_f,line4);
            l_bath_init_pos[j] = std::stod(line1);
            l_bath_init_vel[j] = std::stod(line2);
            r_bath_init_pos[j] = std::stod(line3);
            r_bath_init_vel[j] = std::stod(line4);
        }
    }
    else{
        std::cout << "FILE FAILED TO OPEN. TERMINATING EXECUTION." << std::endl;
        cleanup();
        return 0;
    }

    // Close the files
    left_noise_pos_f.close();
    left_noise_vel_f.close();
    right_noise_pos_f.close();
    right_noise_vel_f.close();

    std::cout << "Noise data successfully loaded" << std::endl;

   Weights_and_Nodes left_wghts_nds {53};
   Weights_and_Nodes right_wghts_nds {53};



    // Initialize the positions and velicities of the resolved system to zero
    std::cout << "Initializing positions and velocities of the resolved system" << std::endl;
    for (int j = 0; j < n_particles; ++j)
    {
        positions[j] = 0;
    }
    first_vel[0] = 0;
    last_vel[0] = 0;
    
    for (int j = 0; j < n_int_particles; ++j)
    {
        int_velocities[j] = 0;
    }

    // Initialize the time domain values for the beta, a_0(t) and b_0(t)
    std::cout << "Initializing bath auxillary functions" << std::endl;

    // Calculate theta normalization
    theta[0] = 1;
    for (long long int j = 1; j < std::min(max_trap_steps, total_steps); ++j)
    {
        theta[j] = theta_f(j*g_dt, left_wghts_nds); // Doesn't matter left or right here
    }

    trunc_theta_norm = normalize_theta( theta, 1.0 );
    std::cout << "The normalization of beta and theta are " << trunc_beta_norm <<
                " and " << trunc_theta_norm << std::endl;

    
    // Initialize theta on the quarter grid
    theta[0] = 1;
    for (long long int j = 1; j < quarter_grid_steps + 1; ++j)
    {
        theta[j] = theta_f(j * g_dt * 0.25, left_wghts_nds); //Doesnt matter left or right here
    }


    a_0[0] = 0;
    b_0[0] = 1;

    // The 4.0 is for the quarter grid values which we determine first
    for (long long int j = 1; j < quarter_grid_steps + 1; ++j)
    {
        a_0[j] = a_0_f(j*g_dt/4.0);
        b_0[j] = b_0_f(j*g_dt/4.0);
    }


    //for (long long int j = 0; j < total_steps; ++j)
    //{
     //   theta[j] = theta_fine[2*j];
    //}
    
    // Initialize the values for the vectors a and b in the bath for j = 1,..., n_bath_particles
    for (long long int j = 0; j < n_bath_particles - 1; ++j)
    {
        a_bath[j] = 0;
        b_bath[j] = 0;
    }


    // Start the VV integrator for the resolved system dynamics
    double l_noise { 0 };
    double r_noise { 0 };
    l_noise = a_0[0]*l_bath_init_vel[0] + b_0[0]*l_bath_init_pos[0];
    r_noise = a_0[0]*r_bath_init_vel[0] + b_0[0]*r_bath_init_pos[0];

    std::ofstream l_noise_f ("left_noise.csv");
    std::ofstream r_noise_f ("right_noise.csv");
    l_noise_f << std::setprecision(std::numeric_limits<double>::max_digits10);
    r_noise_f << std::setprecision(std::numeric_limits<double>::max_digits10);
    
    
    // Calculate the norm of the truncated theta
    //trunc_theta_norm = normalize_theta( theta, 0.25 );

    // For debugging with noise off only
    //quarter_grid_first_vel[0] = 0.01;
    //positions[2] = 0.05;

    /*

    std::cout << "Computing 8 values on the quarter grid with Forward Euler" << std::endl;
    
    double* temp_quarter_grid_int_vel = new double[n_int_particles * quarter_grid_steps]();
    double* temp_quarter_grid_pos = new double[n_particles * quarter_grid_steps]();

    for (int j = 0; j < n_particles; ++j)
    {
        temp_quarter_grid_pos[j] = quarter_grid_pos[j];
    }

    for (int j = 0; j < n_int_particles; ++j)
    {
        temp_quarter_grid_int_vel[j] = quarter_grid_int_vel[j];
    }
    double temp_quarter_grid_first_vel {quarter_grid_first_vel[0]};
    double temp_quarter_grid_last_vel {quarter_grid_last_vel[0]};

    f_system_finer_grids(quarter_grid_pos, quarter_grid_first_vel, quarter_grid_last_vel, quarter_grid_int_vel, 1, theta, l_noise, r_noise, trunc_theta_norm, 0.25);

    for (int j = 0; j < n_particles; ++j)
    {
        quarter_grid_pos[j] = temp_quarter_grid_pos[j]
            + 0.25 * g_dt * quarter_grid_pos[j];
        store_quarter_grid_pos[n_particles + j] = quarter_grid_pos[j];
        store_quarter_grid_pos[j] = temp_quarter_grid_pos[j];
    }
    for (int j = 0; j < n_int_particles; ++j)
    {
        quarter_grid_int_vel[j] = temp_quarter_grid_int_vel[j]
            + 0.25 * g_dt * quarter_grid_int_vel[j];
        store_quarter_grid_int_vel[n_int_particles + j] = quarter_grid_int_vel[j];
        store_quarter_grid_int_vel[j] = temp_quarter_grid_int_vel[j];
    }
    quarter_grid_first_vel[1] = temp_quarter_grid_first_vel + 0.25 * g_dt * quarter_grid_first_vel[0];
    quarter_grid_first_vel[0] = temp_quarter_grid_first_vel; 
    quarter_grid_last_vel[1] = temp_quarter_grid_last_vel + 0.25 * g_dt * quarter_grid_last_vel[0];
    quarter_grid_last_vel[0] = temp_quarter_grid_last_vel;

    
    update_bath_time_zero_quarter( a_0, b_0, a_bath, b_bath, 0);
    l_noise =  F( a_0, b_0, a_bath, b_bath, l_bath_init_pos, l_bath_init_vel, 1);
    r_noise =  F( a_0, b_0, a_bath, b_bath, r_bath_init_pos, r_bath_init_vel, 1);


    for (int current_quarter_step = 1; current_quarter_step < quarter_grid_steps - 1; ++current_quarter_step)
    {
        for (int j = 0; j < n_particles; ++j)
        {
            temp_quarter_grid_pos[j] = quarter_grid_pos[j];
        }

        for (int j = 0; j < n_int_particles; ++j)
        {
            temp_quarter_grid_int_vel[j] = quarter_grid_int_vel[j];
        }

        temp_quarter_grid_first_vel = quarter_grid_first_vel[current_quarter_step];
        temp_quarter_grid_last_vel = quarter_grid_last_vel[current_quarter_step];
    
        f_system_finer_grids(quarter_grid_pos, quarter_grid_first_vel, quarter_grid_last_vel, quarter_grid_int_vel, current_quarter_step + 1, theta, l_noise, r_noise, trunc_theta_norm, 0.25);

        for (int j = 0; j < n_particles; ++j)
        {
            quarter_grid_pos[j] = temp_quarter_grid_pos[j]
                + 0.25 * g_dt * quarter_grid_pos[j];
            store_quarter_grid_pos[n_particles * (current_quarter_step + 1) + j] = quarter_grid_pos[j];
            store_quarter_grid_pos[n_particles * current_quarter_step + j] = temp_quarter_grid_pos[j];
        }
        for (int j = 0; j < n_int_particles; ++j)
        {
            quarter_grid_int_vel[j] = temp_quarter_grid_int_vel[j]
                + 0.25 * g_dt * quarter_grid_int_vel[j];
            store_quarter_grid_int_vel[n_int_particles * (current_quarter_step + 1) + j] = quarter_grid_int_vel[j];
            store_quarter_grid_int_vel[n_int_particles * current_quarter_step + j] = temp_quarter_grid_int_vel[j];
        }
        quarter_grid_first_vel[current_quarter_step + 1] = temp_quarter_grid_first_vel + 0.25 * g_dt * quarter_grid_first_vel[current_quarter_step];
        quarter_grid_first_vel[current_quarter_step] = temp_quarter_grid_first_vel;
        quarter_grid_last_vel[current_quarter_step + 1] = temp_quarter_grid_last_vel + 0.25 * g_dt * quarter_grid_last_vel[current_quarter_step];
        quarter_grid_last_vel[current_quarter_step] = temp_quarter_grid_last_vel;

        

        update_bath_quarter( a_0, b_0, a_bath, b_bath, g_dt * 0.25 * current_quarter_step);
        l_noise =  F( a_0, b_0, a_bath, b_bath, l_bath_init_pos, l_bath_init_vel, current_quarter_step + 1);
        //std::cout << "Check bath order of convergence at time " << g_dt * 0.25 * (current_quarter_step + 1) << " : " << l_noise << std::endl;
        r_noise =  F( a_0, b_0, a_bath, b_bath, r_bath_init_pos, r_bath_init_vel, current_quarter_step + 1);


    }

    delete[] temp_quarter_grid_int_vel;
    temp_quarter_grid_int_vel = nullptr;
    delete[] temp_quarter_grid_pos;
    temp_quarter_grid_pos = nullptr;


    // Re-initialize the noise
    std::cout << "Reinitializing the noise for half grid computation" << std::endl;

        // Double check this
    theta[0] = 1;
    for (long long int j = 1; j < half_grid_steps + 1; ++j)
    {
        theta[j] = theta_f(j * g_dt * 0.5);
    }

    //trunc_theta_norm = normalize_theta( theta, 0.5 );

    a_0[0] = 0;
    b_0[0] = 1;
    for (long long int j = 1; j < half_grid_steps + 1; ++j)
    {
        a_0[j] = a_0_f(j*g_dt/2.0);
        b_0[j] = b_0_f(j*g_dt/2.0);
    }

    //for (long long int j = 0; j < total_steps; ++j)
    //{
     //   theta[j] = theta_fine[2*j];
    //}
    
    // Initialize the values for the vectors a and b in the bath for j = 1,..., n_bath_particles
    for (long long int j = 0; j < n_bath_particles - 1; ++j)
    {
        a_bath[j] = 0;
        b_bath[j] = 0;
    }

    // Start the VV integrator for the resolved system dynamics
    l_noise = a_0[0]*l_bath_init_vel[0] + b_0[0]*l_bath_init_pos[0];
    r_noise = a_0[0]*r_bath_init_vel[0] + b_0[0]*r_bath_init_pos[0];

    std::cout << "Computing the positions and velocities with forward Euler on the half grid" << std::endl;


    double* temp_half_grid_int_vel = new double[n_int_particles * half_grid_steps]();
    double* temp_half_grid_pos = new double[n_particles * half_grid_steps]();

    for (int j = 0; j < n_particles; ++j)
    {
        temp_half_grid_pos[j] = half_grid_pos[j];
    }

    for (int j = 0; j < n_int_particles; ++j)
    {
        temp_half_grid_int_vel[j] = half_grid_int_vel[j];
    }
    double temp_half_grid_first_vel {half_grid_first_vel[0]};
    double temp_half_grid_last_vel {half_grid_last_vel[0]};

    f_system_finer_grids(half_grid_pos, half_grid_first_vel, half_grid_last_vel, half_grid_int_vel, 1, theta, l_noise, r_noise, trunc_theta_norm, 0.5);

    for (int j = 0; j < n_particles; ++j)
    {
        half_grid_pos[j] = temp_half_grid_pos[j]
            + 0.5 * g_dt * half_grid_pos[j];
        store_half_grid_pos[n_particles + j] = half_grid_pos[j];
        store_half_grid_pos[j] = temp_half_grid_pos[j];
    }
    for (int j = 0; j < n_int_particles; ++j)
    {
        half_grid_int_vel[j] = temp_half_grid_int_vel[j]
            + 0.5 * g_dt * half_grid_int_vel[j];
        store_half_grid_int_vel[n_int_particles + j] = half_grid_int_vel[j];
        store_half_grid_int_vel[j] = temp_half_grid_int_vel[j];
    }
    half_grid_first_vel[1] = temp_half_grid_first_vel + 0.5 * g_dt * half_grid_first_vel[0];
    half_grid_first_vel[0] = temp_half_grid_first_vel;
    half_grid_last_vel[1] = temp_half_grid_last_vel + 0.5 * g_dt * half_grid_last_vel[0];
    half_grid_last_vel[0] = temp_half_grid_last_vel;


    update_bath_time_zero( a_0, b_0, a_bath, b_bath, 0);
    l_noise =  F( a_0, b_0, a_bath, b_bath, l_bath_init_pos, l_bath_init_vel, 1);
    r_noise =  F( a_0, b_0, a_bath, b_bath, r_bath_init_pos, r_bath_init_vel, 1);


    for (int current_half_step = 1; current_half_step < half_grid_steps - 1; ++current_half_step)
    {
        for (int j = 0; j < n_particles; ++j)
        {
            temp_half_grid_pos[j] = half_grid_pos[j];
        }

        for (int j = 0; j < n_int_particles; ++j)
        {
            temp_half_grid_int_vel[j] = half_grid_int_vel[j];
        }

        temp_half_grid_first_vel = half_grid_first_vel[current_half_step];
        temp_half_grid_last_vel = half_grid_last_vel[current_half_step];
    
        f_system_finer_grids(half_grid_pos, half_grid_first_vel, half_grid_last_vel, half_grid_int_vel, current_half_step + 1, theta, l_noise, r_noise, trunc_theta_norm, 0.5);

        for (int j = 0; j < n_particles; ++j)
        {
            half_grid_pos[j] = temp_half_grid_pos[j]
                + 0.5 * g_dt * half_grid_pos[j];
            store_half_grid_pos[n_particles * (current_half_step + 1) + j] = half_grid_pos[j];
            store_half_grid_pos[n_particles * current_half_step + j] = temp_half_grid_pos[j];
        }
        for (int j = 0; j < n_int_particles; ++j)
        {
            half_grid_int_vel[j] = temp_half_grid_int_vel[j]
                + 0.5 * g_dt * half_grid_int_vel[j];
            store_half_grid_int_vel[n_int_particles * (current_half_step + 1) + j] = half_grid_int_vel[j];
            store_half_grid_int_vel[n_int_particles * current_half_step + j] = temp_half_grid_int_vel[j];
        }
        half_grid_first_vel[current_half_step + 1] = temp_half_grid_first_vel + 0.5 * g_dt * half_grid_first_vel[current_half_step];
        half_grid_first_vel[current_half_step] =temp_half_grid_first_vel;
        half_grid_last_vel[current_half_step + 1] = temp_half_grid_last_vel + 0.5 * g_dt * half_grid_last_vel[current_half_step];
        half_grid_last_vel[current_half_step]= temp_half_grid_last_vel;

        update_bath( a_0, b_0, a_bath, b_bath, g_dt * 0.5 * current_half_step);
        l_noise =  F( a_0, b_0, a_bath, b_bath, l_bath_init_pos, l_bath_init_vel, current_half_step + 1);
        r_noise =  F( a_0, b_0, a_bath, b_bath, r_bath_init_pos, r_bath_init_vel, current_half_step + 1);


    }

    delete[] temp_half_grid_int_vel;
    temp_half_grid_int_vel = nullptr;
    delete[] temp_half_grid_pos;
    temp_half_grid_pos = nullptr;

    // Richardson Extrapolation on the half grid.

    for (int j = 0; j < half_grid_steps; ++j)
    {
        half_grid_first_vel[j] = (2.0 * quarter_grid_first_vel[2*j] - half_grid_first_vel[j] );
        half_grid_last_vel[j] = (2.0 * quarter_grid_last_vel[2*j] - half_grid_last_vel[j] );
        for (int i = 0; i < n_particles; ++i)
        {
            store_half_grid_pos[n_particles * j + i] = (2.0 * store_quarter_grid_pos[2 * n_particles * j + i] - store_half_grid_pos[n_particles * j + i]);
        }
        for (int i = i; i < n_int_particles; ++i)
        {
            store_half_grid_int_vel[n_int_particles * j + i] = (2.0 * store_quarter_grid_int_vel[2 * n_int_particles * j + i] - store_half_grid_int_vel[n_int_particles * j + i]);
        }
    }

    std::cout << "CHECK VELOCITY FOR RICHARDSON EXTRAPOLATION at time " << g_dt * 0.5 * 1 << " : " << store_half_grid_pos[n_particles * 1 + n_particles - 1] << std::endl;


    // Re-initialize the noise
    std::cout << "Reinitializing the noise for half RK3 computation" << std::endl;

    // Double check this
    theta[0] = 1;
    for (long long int j = 1; j < total_steps; ++j)
    {
        theta[j] = theta_f(j*g_dt);
    }

    trunc_theta_norm = normalize_theta( theta, 1.0 );

    for (long long int j = 0; j < total_steps; ++j)
    {
        shifted_theta[j] = theta_f(g_dt*(j + 0.5));
    }


    a_0[0] = 0;
    b_0[0] = 1;
    for (long long int j = 1; j < 2*total_steps - 1; ++j)
    {
        a_0[j] = a_0_f(j*g_dt/2.0);
        b_0[j] = b_0_f(j*g_dt/2.0);
    }

    //for (long long int j = 0; j < total_steps; ++j)
    //{
     //   theta[j] = theta_fine[2*j];
    //}
    
    // Initialize the values for the vectors a and b in the bath for j = 1,..., n_bath_particles
    for (long long int j = 0; j < n_bath_particles - 1; ++j)
    {
        a_bath[j] = 0;
        b_bath[j] = 0;
    }

    // Start the VV integrator for the resolved system dynamics
    l_noise = a_0[0]*l_bath_init_vel[0] + b_0[0]*l_bath_init_pos[0];
    r_noise = a_0[0]*r_bath_init_vel[0] + b_0[0]*r_bath_init_pos[0];

    for (int j = 0; j < n_particles; ++j)
    {
        positions[j] = store_half_grid_pos[n_particles * (half_grid_steps - 1) + j];
    }
    for (int j = 0; j < n_int_particles; ++j)
    {
        int_velocities[j] = store_half_grid_int_vel[n_int_particles *(half_grid_steps - 1) + j];
    }

    for (int j = 0; j  < half_grid_steps/2 + 1; ++j)
    {
        first_vel[j] = half_grid_first_vel[2*j];
        last_vel[j] = half_grid_last_vel[2*j];
    }

    update_bath_time_zero(a_0, b_0, a_bath, b_bath, 0);
    l_noise =  F( a_0, b_0, a_bath, b_bath, l_bath_init_pos, l_bath_init_vel, 1);
    r_noise =  F( a_0, b_0, a_bath, b_bath, r_bath_init_pos, r_bath_init_vel, 1);

    for (int j = 1; j < half_grid_steps - 1; ++j )
    {
        update_bath(a_0, b_0, a_bath, b_bath, j * 0.5 * g_dt);
        l_noise =  F( a_0, b_0, a_bath, b_bath, l_bath_init_pos, l_bath_init_vel, j + 1);
        r_noise =  F( a_0, b_0, a_bath, b_bath, r_bath_init_pos, r_bath_init_vel, j + 1);
    }
    


    // Loop until equilibriation time
    std::cout << "Beginning the dynamics" << std::endl;
    
    /*
    */
    RK3_system_time_zero(positions,first_vel, last_vel, int_velocities, 0, theta, l_noise, r_noise,
                            trunc_theta_norm, a_0, b_0, a_bath, b_bath, l_bath_init_pos, 
                            l_bath_init_vel, r_bath_init_pos, r_bath_init_vel, shifted_theta,
                            left_wghts_nds, right_wghts_nds);

    
    // Loop until equilibriation time
    std::cout << "Beginning the dynamics after the initial step" << std::endl;\
                        

    //for (long long int current_step = 2; current_step < steps_to_equi; ++current_step)
    for (long long int current_step = 1; current_step < steps_to_equi; ++current_step)
    {

        // Print Statements for debugging
        if (current_step % dt_inv == 0)
        {
            std::cout << "The Current simulation time is " << current_step*g_dt << std::endl;
            std::cout << "The position of the first particle at (matlab) step " << current_step + 1 << " is " << positions[0] << std::endl;
            std::cout << "The velocity of the first particle at (matlab) step " << current_step + 1 << " is " << first_vel[current_step] << std::endl;
            std::cout << "The position of the last particle at (matlab) step " << current_step + 1 << " is " << positions[n_particles -1] << std::endl;
            std::cout << "The velocity of the last particle at (matlab) step " << current_step + 1 << " is " << last_vel[current_step] << std::endl;
            //std::cout << "The last particle has position " << last_pos[current_step] << " and velocity " << velocities[(current_step % 2)][n_particles-1] << std::endl;
        }

        RK3_system(positions,first_vel, last_vel, int_velocities, current_step, theta, l_noise, r_noise,
                            trunc_theta_norm, a_0, b_0, a_bath, b_bath, l_bath_init_pos, 
                            l_bath_init_vel, r_bath_init_pos, r_bath_init_vel, shifted_theta,
                            left_wghts_nds, right_wghts_nds);


    }


    std::cout << "Reached Predefined Equilibriation time" << std::endl;
    long long int n_sample { 0 };
    // Variables for saving local kinetic temperature
    std::string avg_vel_sqrd_string = "avg_vel_sqrd_";
    std::string csv = "_.csv";
    int loc_temp_file_numb {0};

    for (long long int current_step = steps_to_equi; current_step < total_steps -1; ++current_step)
    {

        // Print Statements for debugging
        if (current_step % dt_inv == 0)
        {
            std::cout << "The Current simulation time is " << current_step*g_dt << std::endl;
            std::cout << "The position of the first particle at (matlab) step " << current_step + 1 << " is " << positions[0] << std::endl;
            std::cout << "The velocity of the first particle at (matlab) step " << current_step + 1 << " is " << first_vel[current_step] << std::endl;
            std::cout << "The position of the last particle at (matlab) step " << current_step + 1 << " is " << positions[n_particles -1] << std::endl;
            std::cout << "The velocity of the last particle at (matlab) step " << current_step + 1 << " is " << last_vel[current_step] << std::endl;
            //std::cout << "The last particle has position " << last_pos[current_step] << " and velocity " << velocities[(current_step % 2)][n_particles-1] << std::endl;
        }


        // Sample the position and velocities of the particles
        // Sample once every 10 steps
        if (current_step % sample_freq == 0)
        {
            vel_samples[n_sample*n_particles] = first_vel[current_step];
            pos_samples[n_sample*n_particles] = positions[0];
            for (int j = 1; j < n_particles-1; ++j)
            {
                vel_samples[n_sample*n_particles + j] = int_velocities[j - 1];
                pos_samples[n_sample*n_particles + j] = positions[j];
            }
            vel_samples[n_sample*n_particles + n_particles -1] = last_vel[current_step];
            pos_samples[n_sample*n_particles + n_particles -1] = positions[n_particles -1];
            n_sample += 1;
        }


        // Save updated local temperature every 1000 seconds after equilibriation
            if (current_step % (dt_inv * 1000) == 0)
            {
                std::cout << "The current time is " << current_step * g_dt << std::endl;
                std::cout << "Calculating local kinetic energy" << std::endl;
                double avg_vel_sqrd[n_particles] {}; // Initialized to all zeros.
                for (int long long i = 0; i < (n_sample-1); ++i)
                {
                    for (int j = 0; j < n_particles; ++j)
                    {
                        avg_vel_sqrd[j] += vel_samples[i*n_particles + j]*vel_samples[i*n_particles + j];
                    }
                }
                for ( int i = 0; i < n_particles; ++i)
                {
                    avg_vel_sqrd[i] = avg_vel_sqrd[i]/(n_sample-1);
                }


                // Saves the local kinetic energy
                std::ofstream avg_vel_sqrd_f ( avg_vel_sqrd_string + std::to_string(loc_temp_file_numb) + ".csv");
                avg_vel_sqrd_f << std::setprecision(std::numeric_limits<double>::max_digits10);
                if (avg_vel_sqrd_f.is_open())
                {
                    for (const double &element : avg_vel_sqrd)
                    {
                        avg_vel_sqrd_f << element << ",";
                    }
                    avg_vel_sqrd_f << std::endl;
                }
                else
                {
                    std::cout<< "ERROR: FAILED TO WRITE LOCAL KINETIC ENERGY";
                    cleanup();
                    return 0;
                }
                avg_vel_sqrd_f.close();

                loc_temp_file_numb += 1;
            }

        RK3_system(positions,first_vel, last_vel, int_velocities, current_step, theta, l_noise, r_noise,
                trunc_theta_norm, a_0, b_0, a_bath, b_bath, l_bath_init_pos, 
                l_bath_init_vel, r_bath_init_pos, r_bath_init_vel, shifted_theta,
                left_wghts_nds, right_wghts_nds);

    }


    std::cout<< "The data is currently at time " << g_dt*(total_steps -1) << std::endl; 
    
    std::cout << "Dynamics complete" <<std::endl;
    // Compute the local kinetic energy as a local temperature substitute

    l_noise_f.close();
    r_noise_f.close();



    if ((total_steps - 1) % sample_freq == 0)
        {
            vel_samples[n_sample*n_particles] = first_vel[(total_steps - 1)];
            pos_samples[n_sample*n_particles] = positions[0];
            for (int j = 1; j < n_particles-1; ++j)
            {
                vel_samples[n_sample*n_particles + j] = int_velocities[j - 1];
                pos_samples[n_sample*n_particles + j] = positions[j];
            }
            vel_samples[n_sample*n_particles + n_particles -1] = last_vel[(total_steps - 1)];
            pos_samples[n_sample*n_particles + n_particles -1] = positions[n_particles -1];
            n_sample += 1;
        }


    std::cout << "Calculating local kinetic energy" << std::endl;
    double avg_vel_sqrd[n_particles] {}; // Initialized to all zeros.
    for (int long long i = 0; i < total_samples; ++i)
    {
        for (int j = 0; j < n_particles; ++j)
        {
            avg_vel_sqrd[j] += vel_samples[i*n_particles + j]*vel_samples[i*n_particles + j];
        }
    }
    for ( int i = 0; i < n_particles; ++i)
    {
        avg_vel_sqrd[i] = avg_vel_sqrd[i]/total_samples;
    }

    std::cout << "Writing files to disk" << std::endl;
    // Saves the Samples to a file to be imported in MATLAB
    std::ofstream vel_samp ("vel_samp.csv");
    std::ofstream pos_samp ("pos_samp.csv");
    vel_samp << std::setprecision(std::numeric_limits<double>::max_digits10);
    pos_samp << std::setprecision(std::numeric_limits<double>::max_digits10);

    if (vel_samp.is_open() && pos_samp.is_open())
    {
        for ( long long int i = 0; i < total_samples; ++i)
        {
            for (int j = 0; j < n_particles; ++j)
            {
                vel_samp << vel_samples[i*n_particles + j] << ",";
                pos_samp << pos_samples[i*n_particles + j] << ",";
            }
            vel_samp << std::endl;
            pos_samp << std::endl;
        }
    }
    else
    {
        std::cout << "ERROR: FAILED TO WRITE SAMPLE FILES!" << std::endl;
        cleanup();
        return 0;
    }
    vel_samp.close();
    pos_samp.close();

    // Saves the local kinetic energy
    std::ofstream avg_vel_sqrd_f ("avg_vel_sqrd_final.csv");
    avg_vel_sqrd_f << std::setprecision(std::numeric_limits<double>::max_digits10);
    if (avg_vel_sqrd_f.is_open())
    {
        for (const double &element : avg_vel_sqrd)
        {
            avg_vel_sqrd_f << element << ",";
        }
        avg_vel_sqrd_f << std::endl;
    }
    else
    {
        std::cout<< "ERROR: FAILED TO WRITE LOCAL KINETIC ENERGY";
        cleanup();
        return 0;
    }
    avg_vel_sqrd_f.close();

    std::cout << "Execution Sucessfull. Cleaning up" <<std::endl;
    cleanup();
    return 0;
}

