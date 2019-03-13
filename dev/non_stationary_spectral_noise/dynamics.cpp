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

    // Load the noise for the baths
    std::ifstream left_noise_f;
    std::ifstream right_noise_f;
    left_noise_f.open("left_noise.dat", std::ios::binary | std::ios::ate);
    right_noise_f.open("right_noise.dat", std::ios::binary | std::ios::ate);

    double* left_noise;
    double* right_noise;

    if ( left_noise_f.is_open() && right_noise_f.is_open() )
    {
        std::streampos left_noise_size { left_noise_f.tellg() };
        std::streampos right_noise_size { right_noise_f.tellg() };

        char* memblock_left_noise = new char[left_noise_size];
        char* memblock_right_noise = new char[right_noise_size];

        left_noise_f.seekg(0, std::ios::beg);
        right_noise_f.seekg(0, std::ios::beg);

        left_noise_f.read(memblock_left_noise, left_noise_size );
        right_noise_f.read(memblock_right_noise, right_noise_size);

        left_noise_f.close();
        right_noise_f.close();

        // Reinterpret as doubles
        left_noise = (double*) memblock_left_noise;
        right_noise = (double*) memblock_right_noise;

        std::cout << "Loaded bath noise " << std::endl;
    }
    else
    {
        std::cerr << "Failed to open noise data files" << std::endl;
    }
    


   Weights_and_Nodes left_wghts_nds {55};
   Weights_and_Nodes right_wghts_nds {55};



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
    

    initialize_positions(positions, n_particles);
    initialize_velocities(int_velocities, first_vel, last_vel, n_int_particles);
    first_initial_pos = positions[0];
    last_initial_pos = positions[n_particles -1];

    /*
    for (int jj = 0; jj < n_particles; ++jj)
    {
        std::cout << "Positions: " << positions[jj]  << std::endl;
    }
    std::cout << "First initial Pos " << first_initial_pos << std::endl;
    std::cout << "Last Initial pos " << last_initial_pos << std::endl;
    std::cout << "First Vel " << first_vel[0] << std::endl;
    for (int jj = 0; jj < n_int_particles; ++jj)
    {
        std::cout << "Velocities " << int_velocities[jj] << std::endl;
    }
    std::cout << "Last Vel " << last_vel[0] << std::endl;
    */
    compute_mean(positions, n_particles);
    compute_variance(positions, n_particles);
    compute_mean(int_velocities, n_int_particles);
    compute_variance(int_velocities, n_int_particles);
    //return 0;



    // Initialize the time domain values for the beta, a_0(t) and b_0(t)
    std::cout << "Initializing bath auxillary functions" << std::endl;

    // Calculate theta normalization
    theta[0] = 1;
    for (long long int j = 1; j < std::min(max_trap_steps, total_steps); ++j)
    {
        theta[j] = theta_f(j*g_dt, left_wghts_nds); // Doesn't matter left or right here
    }

    // trunc_theta_norm = normalize_theta( theta, 1.0 );
    trunc_theta_norm = 1.0; // Theta isnt truncated in the strong solver
    std::cout << "The normalization of beta and theta are " << trunc_beta_norm <<
                " and " << trunc_theta_norm << std::endl;

    
    // Initialize theta on the quarter grid
    theta[0] = 1;
    for (long long int j = 1; j < quarter_grid_steps + 1; ++j)
    {
        theta[j] = theta_f(j * g_dt * 0.25, left_wghts_nds); //Doesnt matter left or right here
    }


    // Loop until equilibriation time
    std::cout << "Beginning the dynamics" << std::endl;
    
    /*
    */
    RK3_system_time_zero(positions,first_vel, last_vel, int_velocities, 0, theta, left_noise, right_noise,
                            trunc_theta_norm, shifted_theta,
                            left_wghts_nds, right_wghts_nds, first_initial_pos, last_initial_pos);

    
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

        RK3_system(positions,first_vel, last_vel, int_velocities, current_step, theta, left_noise, right_noise,
                            trunc_theta_norm, shifted_theta,
                            left_wghts_nds, right_wghts_nds, first_initial_pos, last_initial_pos);


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

        RK3_system(positions,first_vel, last_vel, int_velocities, current_step, theta, left_noise, right_noise,
                trunc_theta_norm, shifted_theta,
                left_wghts_nds, right_wghts_nds, first_initial_pos, last_initial_pos);

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

    delete[] left_noise;
    left_noise = nullptr;
    delete[] right_noise;
    right_noise = nullptr;

    return 0;
}

