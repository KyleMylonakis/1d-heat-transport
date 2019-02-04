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


    std::cout << std::setprecision(std::numeric_limits<double>::max_digits10);
    
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

    std::cout << "Max steps has value " << total_steps << std::endl;
    std::cout << "Max trap steps has value " << max_trap_steps << std::endl;
    

    // Initialize the positions and velicities of the resolved system to zero
    std::cout << "Initializing positions and velocities of the resolved system" << std::endl;
    /*
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
    */
    
    initialize_positions(positions, n_particles);
    initialize_velocities(int_velocities, first_vel, last_vel, n_int_particles);
    first_initial_pos = positions[0];
    last_initial_pos = positions[n_particles -1];
    compute_mean(positions, n_particles);
    compute_variance(positions, n_particles);
    compute_mean(int_velocities, n_int_particles);
    compute_variance(int_velocities, n_int_particles);
    //return 0;
    
    
    // Generate and save noise for testing purposes
    /*
    Noise left_noise {0.2, t_max, g_dt};
    std::ofstream test_noise_f ("test_noise.csv");
    test_noise_f << std::setprecision(std::numeric_limits<double>::max_digits10);
    if (test_noise_f.is_open())
    {

        for (long long int i = 0; i < left_noise.length(); ++i)
        {
            test_noise_f << left_noise[i] << ",";
        }
        test_noise_f << std::endl;
    }
    else
    {
        std::cout<< "ERROR: FAILED TO WRITE TEST NOISE CSV";
        cleanup();
        return 0;
    }
    test_noise_f.close();
    */


   Weights_and_Nodes left_wghts_nds {55};
   Weights_and_Nodes right_wghts_nds {55};


    // Calculate theta normalization
    theta[0] = 1;
    for (long long int j = 1; j < total_steps; ++j)
    {
        // It doesn't matter which eight or node i use here
        theta[j] = theta_f(j*g_dt, left_wghts_nds);
    }

    trunc_theta_norm = normalize_theta( theta, 1.0 );
    
    std::cout << "The normalization of theta is " << trunc_theta_norm << std::endl;
    std::cout << "The integral is truncated at " << Tc << " = " << (max_trap_steps)*g_dt << std::endl;
    std::cout << "Max trap steps is " << max_trap_steps << std::endl;
    //return 0;

    
    std::cout << "Constructing Noise Objects on the eigth grid" << std::endl;

    std::cout << "Constructing Left Noise" << std::endl;
    Noise left_eigth_noise { (left_temp/trunc_theta_norm )/g_k_sqrd ,t_max,g_dt/8.0, left_wghts_nds};
    //Noise left_eigth_noise { (left_temp/trunc_theta_norm ) ,t_max,g_dt/8.0, left_wghts_nds};


    std::cout << "Constructing Noise for richardson extrapolation" << std::endl;
    for (int i = 0; i < eigth_grid_steps; ++i)
    {
        eigth_grid_left_noise[i] = left_eigth_noise[i];
    }
    for (int i = 0; i < quarter_grid_steps; ++i)
    {
        quarter_grid_left_noise[i] = eigth_grid_left_noise[2*i];
    }
    for (int i =0; i < half_grid_steps; ++i)
    {
        half_grid_left_noise[i] = quarter_grid_left_noise[2*i];
    }

    std::cout << "Initializing the noise on the full half grid" << std::endl;
    for (long long int i = 0; i < 2*total_steps - 1; ++i)
    {
        left_noise_half_grid_all[i] = left_eigth_noise[4*i];
    }

    left_eigth_noise.~Noise();

    std::cout << "Constructing Right noise" << std::endl;
    Noise right_eigth_noise  { (right_temp/trunc_theta_norm)/g_k_sqrd, t_max, g_dt/8.0, right_wghts_nds};
    //Noise right_eigth_noise  { (right_temp/trunc_theta_norm), t_max, g_dt/8.0, right_wghts_nds};


    std::cout << "Constructing Noise for richardson extrapolation" << std::endl;
    for (int i = 0; i < eigth_grid_steps; ++i)
    {
        eigth_grid_right_noise[i] = right_eigth_noise[i];
    }
    for (int i = 0; i < quarter_grid_steps; ++i)
    {
        quarter_grid_right_noise[i] = eigth_grid_right_noise[2*i];
    }
    for (int i =0; i < half_grid_steps; ++i)
    {
        half_grid_right_noise[i] = quarter_grid_right_noise[2*i];
    }

    std::cout << "Initializing the noise on the full half grid" << std::endl;
    for (long long int i = 0; i < 2*total_steps - 1; ++i)
    {
        right_noise_half_grid_all[i] = right_eigth_noise[4*i];
    }

    right_eigth_noise.~Noise();

    /*

    /*
    std::cout << "SETTING DEBUG VALUES FOR RICHARDSON EXTRAPOLATION" << std::endl;

    for (int i = 0; i < eigth_grid_steps; ++i)
    {
        eigth_grid_left_noise[i] = 0;
        eigth_grid_right_noise[i] = 0;
    }
    for (int i = 0; i < quarter_grid_steps; ++i)
    {
        quarter_grid_left_noise[i] = 0;
        quarter_grid_right_noise[i] = 0;
    }
    for (int i =0; i < half_grid_steps; ++i)
    {
        half_grid_left_noise[i] = 0;
        half_grid_right_noise[i] = 0;
    }

    double DEBUG_VEL {0.1};
    quarter_grid_first_vel[0] = DEBUG_VEL;
    half_grid_first_vel[0] = DEBUG_VEL;
    eigth_grid_first_vel[0] = DEBUG_VEL;
    for (int j = 0; j < n_int_particles; ++j)
    {
        eigth_grid_int_vel[j] = DEBUG_VEL;
        quarter_grid_int_vel[j] = DEBUG_VEL;
        half_grid_int_vel[j] = DEBUG_VEL;
    }
    quarter_grid_last_vel[0] = DEBUG_VEL;
    half_grid_last_vel[0] = DEBUG_VEL;
    eigth_grid_last_vel[0] = DEBUG_VEL;
    */

   /*
    std::cout << "Computing 16 values on the eigth grid with Forward Euler" << std::endl;
    
    theta[0] = 1;
    for (long long int j = 1; j < eigth_grid_steps + 1; ++j)
    {
        theta[j] = theta_f(j * g_dt * 0.125);
    }

   
    double* temp_eigth_grid_int_vel = new double[n_int_particles * eigth_grid_steps]();
    double* temp_eigth_grid_pos = new double[n_particles * eigth_grid_steps]();

    for (int j = 0; j < n_particles; ++j)
    {
        temp_eigth_grid_pos[j] = eigth_grid_pos[j];
    }

    for (int j = 0; j < n_int_particles; ++j)
    {
        temp_eigth_grid_int_vel[j] = eigth_grid_int_vel[j];
    }
    double temp_eigth_grid_first_vel {eigth_grid_first_vel[0]};
    double temp_eigth_grid_last_vel {eigth_grid_last_vel[0]};

    
    f_system_finer_grids(eigth_grid_pos, eigth_grid_first_vel, eigth_grid_last_vel, eigth_grid_int_vel, 1, theta, eigth_grid_left_noise, eigth_grid_right_noise, trunc_theta_norm, 0.125);

    for (int j = 0; j < n_particles; ++j)
    {
        eigth_grid_pos[j] = temp_eigth_grid_pos[j]
            + 0.125 * g_dt * eigth_grid_pos[j];
        store_eigth_grid_pos[n_particles + j] = eigth_grid_pos[j];
        store_eigth_grid_pos[j] = temp_eigth_grid_pos[j];
    }
    for (int j = 0; j < n_int_particles; ++j)
    {
        eigth_grid_int_vel[j] = temp_eigth_grid_int_vel[j]
            + 0.125 * g_dt * eigth_grid_int_vel[j];
        store_eigth_grid_int_vel[n_int_particles + j] = eigth_grid_int_vel[j];
        store_eigth_grid_int_vel[j] = temp_eigth_grid_int_vel[j];
    }
    eigth_grid_first_vel[1] = temp_eigth_grid_first_vel + 0.125 * g_dt * eigth_grid_first_vel[0];
    eigth_grid_first_vel[0] = temp_eigth_grid_first_vel;
    eigth_grid_last_vel[1] = temp_eigth_grid_last_vel + 0.125 * g_dt * eigth_grid_last_vel[0];
    eigth_grid_last_vel[0] = temp_eigth_grid_last_vel;

    for (int current_eigth_step = 1; current_eigth_step < eigth_grid_steps - 1; ++current_eigth_step)
    {
        for (int j = 0; j < n_particles; ++j)
        {
            temp_eigth_grid_pos[j] = eigth_grid_pos[j];
        }

        for (int j = 0; j < n_int_particles; ++j)
        {
            temp_eigth_grid_int_vel[j] = eigth_grid_int_vel[j];
        }

        temp_eigth_grid_first_vel = eigth_grid_first_vel[current_eigth_step];
        temp_eigth_grid_last_vel = eigth_grid_last_vel[current_eigth_step];
    
        f_system_finer_grids(eigth_grid_pos, eigth_grid_first_vel, eigth_grid_last_vel, eigth_grid_int_vel, current_eigth_step + 1, theta, eigth_grid_left_noise, eigth_grid_right_noise, trunc_theta_norm, 0.125);

        for (int j = 0; j < n_particles; ++j)
        {
            eigth_grid_pos[j] = temp_eigth_grid_pos[j]
                + 0.125 * g_dt * eigth_grid_pos[j];
            store_eigth_grid_pos[n_particles * (current_eigth_step + 1) + j] = eigth_grid_pos[j];
            store_eigth_grid_pos[n_particles * current_eigth_step + j] = temp_eigth_grid_pos[j];
        }
        for (int j = 0; j < n_int_particles; ++j)
        {
            eigth_grid_int_vel[j] = temp_eigth_grid_int_vel[j]
                + 0.125 * g_dt * eigth_grid_int_vel[j];
            store_eigth_grid_int_vel[n_int_particles * (current_eigth_step + 1) + j] = eigth_grid_int_vel[j];
            store_eigth_grid_int_vel[n_int_particles * current_eigth_step + j] = temp_eigth_grid_int_vel[j];
        }
        eigth_grid_first_vel[current_eigth_step + 1] = temp_eigth_grid_first_vel + 0.125 * g_dt * eigth_grid_first_vel[current_eigth_step];
        eigth_grid_first_vel[current_eigth_step] =temp_eigth_grid_first_vel;
        eigth_grid_last_vel[current_eigth_step + 1] = temp_eigth_grid_last_vel + 0.125 * g_dt * eigth_grid_last_vel[current_eigth_step];
        eigth_grid_last_vel[current_eigth_step]= temp_eigth_grid_last_vel;

    }

    delete[] temp_eigth_grid_int_vel;
    temp_eigth_grid_int_vel = nullptr;
    delete[] temp_eigth_grid_pos;
    temp_eigth_grid_pos = nullptr;


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

    f_system_finer_grids(quarter_grid_pos, quarter_grid_first_vel, quarter_grid_last_vel, quarter_grid_int_vel, 1, theta, quarter_grid_left_noise, quarter_grid_right_noise, trunc_theta_norm, 0.25);

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
    
        f_system_finer_grids(quarter_grid_pos, quarter_grid_first_vel, quarter_grid_last_vel, quarter_grid_int_vel, current_quarter_step + 1, theta, quarter_grid_left_noise, quarter_grid_right_noise, trunc_theta_norm, 0.25);

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

    }

    delete[] temp_quarter_grid_int_vel;
    temp_quarter_grid_int_vel = nullptr;
    delete[] temp_quarter_grid_pos;
    temp_quarter_grid_pos = nullptr;


        // Double check this
    theta[0] = 1;
    for (long long int j = 1; j < half_grid_steps + 1; ++j)
    {
        theta[j] = theta_f(j * g_dt * 0.5);
    }
    

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

    f_system_finer_grids(half_grid_pos, half_grid_first_vel, half_grid_last_vel, half_grid_int_vel, 1, theta, half_grid_left_noise, half_grid_right_noise, trunc_theta_norm, 0.5);

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
    
        f_system_finer_grids(half_grid_pos, half_grid_first_vel, half_grid_last_vel, half_grid_int_vel, current_half_step + 1, theta, half_grid_left_noise, half_grid_right_noise, trunc_theta_norm, 0.5);

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

    }

    delete[] temp_half_grid_int_vel;
    temp_half_grid_int_vel = nullptr;
    delete[] temp_half_grid_pos;
    temp_half_grid_pos = nullptr;

    // Richardson Extrapolation on the half grid.

    double* ptr_temp_quarter_grid_first_vel = new double[quarter_grid_steps]();
    double* ptr_temp_quarter_grid_last_vel = new double[quarter_grid_steps]();
    double* temp_store_quarter_grid_pos = new double[n_particles*quarter_grid_steps]();
    double* temp_store_quarter_grid_int_vel = new double[n_int_particles*quarter_grid_steps]();

    for (int j = 0; j < quarter_grid_steps; ++j)
    {
        ptr_temp_quarter_grid_first_vel[j] = (2.0 * eigth_grid_first_vel[2*j] - quarter_grid_first_vel[j] );
        ptr_temp_quarter_grid_last_vel[j] = (2.0 * eigth_grid_last_vel[2*j] - quarter_grid_last_vel[j] );
        for (int i = 0; i < n_particles; ++i)
        {
            temp_store_quarter_grid_pos[n_particles * j + i] = (2.0 * store_eigth_grid_pos[2 * n_particles * j + i] - store_quarter_grid_pos[n_particles * j + i]);
        }
        for (int i = 0; i < n_int_particles; ++i)
        {
            temp_store_quarter_grid_int_vel[n_int_particles * j + i] = (2.0 * store_eigth_grid_int_vel[2 * n_int_particles * j + i] - store_quarter_grid_int_vel[n_int_particles * j + i]);
        }
    }
    
    std::cout << "CHECK VELOCITY FOR RICHARDSON EXTRAPOLATION QUARTER GRID at time " << g_dt * 0.25 * ((quarter_grid_steps - 1)/mod_time_step) << " : " << ptr_temp_quarter_grid_first_vel[(quarter_grid_steps - 1)/mod_time_step] << std::endl;

    double* ptr_temp_half_grid_first_vel = new double[half_grid_steps]();
    double* ptr_temp_half_grid_last_vel = new double[half_grid_steps]();
    double* temp_store_half_grid_pos = new double[n_particles*half_grid_steps]();
    double* temp_store_half_grid_int_vel = new double[n_int_particles*half_grid_steps]();


    
    for (int j = 0; j < half_grid_steps; ++j)
    {
        ptr_temp_half_grid_first_vel[j] = (2.0 * quarter_grid_first_vel[2*j] - half_grid_first_vel[j] );
        ptr_temp_half_grid_last_vel[j] = (2.0 * quarter_grid_last_vel[2*j] - half_grid_last_vel[j] );
        for (int i = 0; i < n_particles; ++i)
        {
            temp_store_half_grid_pos[n_particles * j + i] = (2.0 * store_quarter_grid_pos[2 * n_particles * j + i] - store_half_grid_pos[n_particles * j + i]);
        }
        for (int i = 0; i < n_int_particles; ++i)
        {
            temp_store_half_grid_int_vel[n_int_particles * j + i] = (2.0 * store_quarter_grid_int_vel[2 * n_int_particles * j + i] - store_half_grid_int_vel[n_int_particles * j + i]);
        }
    }

    std::cout << "CHECK VELOCITY FOR RICHARDSON EXTRAPOLATION HALF GRID at time " << g_dt * 0.5 * ((half_grid_steps - 1)/mod_time_step) << " : " << ptr_temp_half_grid_first_vel[(half_grid_steps - 1)/mod_time_step] << std::endl;

    

    // Get third order from another extrapolation
    
    for (int j = 0; j < half_grid_steps; ++j)
    {
        //std::cout << "The values of the fine and corase grid are :" << ptr_temp_quarter_grid_first_vel[2*j] << " " << ptr_temp_half_grid_first_vel[j] << std::endl;
        //std::cout << "The value of the extrapolation is " << (16.0 * ptr_temp_quarter_grid_first_vel[2*j] - ptr_temp_half_grid_first_vel[j] )/15.0 << std::endl;
        //std::cout << "The value of the extrapolation is " << (8 * ptr_temp_quarter_grid_first_vel[2*j] - ptr_temp_half_grid_first_vel[j] )/7.0 << std::endl;
        /*for (int k = 0; k < half_grid_steps; ++k)
        {
            std::cout << ptr_temp_half_grid_first_vel[k] << std::endl;
        }*/

        /*
        half_grid_first_vel[j] = (4.0 * ptr_temp_quarter_grid_first_vel[2*j] - 1.0*ptr_temp_half_grid_first_vel[j] )/3.0;
        //half_grid_first_vel[j] = (-1.0 * ptr_temp_quarter_grid_first_vel[2*j] + 4.0*ptr_temp_half_grid_first_vel[j] )/3.0;
        half_grid_last_vel[j] = (4.0 * ptr_temp_quarter_grid_last_vel[2*j] - ptr_temp_half_grid_last_vel[j] )/3.0;
        for (int i = 0; i < n_particles; ++i)
        {
            store_half_grid_pos[n_particles * j + i] = (4.0 * temp_store_quarter_grid_pos[2 * n_particles * j + i] - temp_store_half_grid_pos[n_particles * j + i])/3.0;
        }
        for (int i = 0; i < n_int_particles; ++i)
        {
            store_half_grid_int_vel[n_int_particles * j + i] = (4.0 * temp_store_quarter_grid_int_vel[2 * n_int_particles * j + i] - temp_store_half_grid_int_vel[n_int_particles * j + i])/3.0;
        }
    }
    
    std::cout << "CHECK VELOCITY FOR THIRD ORDER RICHARDSON EXTRAPOLATION at time " << g_dt * 0.5 * ((half_grid_steps - 1)/mod_time_step) << " : " << half_grid_first_vel[(half_grid_steps - 1)/mod_time_step] << std::endl;


    //return 0;


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
    //for (long long int j = 0; j < total_steps; ++j)
    //{
     //   theta[j] = theta_fine[2*j];
    //}

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


    delete[] ptr_temp_quarter_grid_first_vel;
    delete[] ptr_temp_quarter_grid_last_vel;
    delete[] temp_store_quarter_grid_pos;
    delete[] temp_store_quarter_grid_int_vel;

    delete[] ptr_temp_half_grid_first_vel;
    delete[] ptr_temp_half_grid_last_vel;
    delete[] temp_store_half_grid_pos;
    delete[] temp_store_half_grid_int_vel;



    ptr_temp_quarter_grid_first_vel = nullptr;
    ptr_temp_quarter_grid_last_vel = nullptr;
    temp_store_quarter_grid_pos = nullptr;
    temp_store_quarter_grid_int_vel = nullptr;


    ptr_temp_half_grid_first_vel = nullptr;
    ptr_temp_half_grid_last_vel = nullptr;
    temp_store_half_grid_pos = nullptr;
    temp_store_half_grid_int_vel = nullptr;


    /*
    // DEBUG VALUES FOR NOISE and SYSTEM
    std::cout << "Initializing the noise on the full half grid" << std::endl;
    for (long long int i = 0; i < 2*total_steps - 1; ++i)
    {
        left_noise_half_grid_all[i] = 0;
        right_noise_half_grid_all[i] = 0;
    }
        for (int j = 0; j < n_particles; ++j)
    {
        positions[j] = 0;
    }
    for (long long int i = 0; i < total_steps; ++i)
    {
        first_vel[i] = 0;
        last_vel[i] = 0;
    }

    first_vel[8/mod_time_step] = 1;
    
    
    for (int j = 0; j < n_int_particles; ++j)
    {
        int_velocities[j] = 0;
    }
    */

    // Loop until equilibriation time
    std::cout << "Beginning the dynamics after the initial Euler steps" << std::endl;
                        

    //for (long long int current_step = 2; current_step < steps_to_equi; ++current_step)
    for (long long int current_step = 0; current_step < steps_to_equi; ++current_step)
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

        RK3_system(positions,first_vel, last_vel, int_velocities, current_step, theta, left_noise_half_grid_all, right_noise_half_grid_all,
                            trunc_theta_norm, shifted_theta, left_wghts_nds,  right_wghts_nds, first_initial_pos, last_initial_pos);


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

        RK3_system(positions,first_vel, last_vel, int_velocities, current_step, theta, left_noise_half_grid_all, right_noise_half_grid_all,
                trunc_theta_norm, shifted_theta, left_wghts_nds, right_wghts_nds, first_initial_pos, last_initial_pos);

    }


    std::cout<< "The data is currently at time " << g_dt*(total_steps -1) << std::endl; 
    
    std::cout << "Dynamics complete" <<std::endl;
    // Compute the local kinetic energy as a local temperature substitute



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

