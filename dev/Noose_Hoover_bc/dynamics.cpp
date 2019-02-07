/*
    * 
    * The main function will simulate a nearest neighbor interacting system
    * sampling at regular intervals after an equilibriation time equal to the 
    * half the final time of the simulation.
    * 
    * Inputs:
    * Int inputs always come before double inputs
    * The main loop will have the following necessary command line parameters:
    * 1) The final time of the simulation
    * 2) The total number of particles
    * 3) Temperature of left bath
    * 4) Noose Hoover Coupling constant for left bath
    * 5) Temperature of right bath
    * 6) Noose Hoover Coupling constant for right bath
    * 7) The time step of the simulation
    * 8) Fixed or free boundary condition
    * 
    * Outputs:
    * Integer representing the exit status of the system
    * 
    * File outputs:
    * Binary file containing the number of particles of the simulation and the number
    * of samples
    * Binary files containing regularly spaced samples of the position and velocity
    * of the system.
*/

// System Includes
#include <iostream>
#include <iomanip>
#include <limits>
#include <chrono>
#include <fstream>
#include <string>

// Custom Includes
#include "parse_input.h"
#include "Particle_system.h"


int main(int argc, char* argv[])
{
    // Set maximum output precision
    std::cout << std::setprecision(std::numeric_limits<double>::max_digits10);

    // Require the number of command line args to be exact
    int num_int_inputs {2};
    int num_double_inputs {5};
    int num_string_inputs {1};
    int num_required_inputs {num_int_inputs + num_double_inputs + num_string_inputs};

    if (num_required_inputs != argc -1 )
    {
        std::cerr << "Error: Not enough input arguments" << std::endl;
        return -1;
    }

    // Parse the inputs
    int t_final {0};
    int n_particles {0};
    double left_temp {0.0};
    double left_NH_coupling_const {0.0};
    double right_temp {0.0};
    double right_NH_coupling_const {0.0};
    double dt {0.0};
    std::string boundary_condition;

    parse_inputs(argv[1], t_final);
    parse_inputs(argv[2], n_particles);
    parse_inputs(argv[3], left_temp);
    parse_inputs(argv[4], left_NH_coupling_const);
    parse_inputs(argv[5], right_temp);
    parse_inputs(argv[6], right_NH_coupling_const);
    parse_inputs(argv[7], dt);
    parse_inputs(argv[8], boundary_condition);

    // Initialize the Particle system
    Particle_system particles { n_particles, 
        left_temp, 
        left_NH_coupling_const, 
        right_temp, 
        right_NH_coupling_const,
        boundary_condition };

    // Initialize the variable controlling the length
    // of the simulation
    int dt_inv {static_cast<int>(1/dt)};
    std::cout << "dt_inv " << dt_inv << std::endl;
    long int total_steps {static_cast<long int>(t_final) * dt_inv };
    long int steps_before_sampling { total_steps /2 };

    // Controls how often logs are reported to standard out
    int report_every {10};

    // Initialize the sample arrays
    double sample_dt {0.1}; // Simulation time between samples
    int sample_every { static_cast<int>(dt_inv*sample_dt) };
    int num_sample_times { (total_steps - steps_before_sampling)/sample_every + 1};

    double* positions_samples = new double[num_sample_times * n_particles]();
    double* velocity_samples = new double[num_sample_times * n_particles]();

    // Begin a timer and start the simulation: simulate without
    // sampling until a preset equilibriation time
    auto start_time {std::chrono::system_clock::now()};

    for (long int ii = 0; ii < steps_before_sampling; ++ii)
    {
        // Report to standard output regularly
        if (ii % (dt_inv*report_every) == 0)
        {
            std::cout << "The system is currently at time " << ii * dt << std::endl;
            std::cout << "The first particle has velocity " << particles.get_velocity(0) << std::endl;
        }
        
        // Advance the dynamics
        particles.RK3(dt);
    }

    // Continue the dynamics, now also sampling at regular intervals
    int samples_written {0};
    for (long int ii = steps_before_sampling; ii < total_steps; ++ii)
    {
        // Report to standard output regualrly
        if (ii % (dt_inv*report_every) == 0)
        {
            std::cout << "The system is currently at time " << ii * dt << std::endl;
            std::cout << "The first particle has velocity " << particles.get_velocity(0) << std::endl;
        }

        // Sample the dynamics regularly
        if ( ii % sample_every == 0 )
        {
            for (int jj = 0; jj < n_particles; ++jj)
            {
                positions_samples[samples_written*n_particles +jj ] = particles.get_position(jj);
                velocity_samples[ samples_written*n_particles +jj ] = particles.get_velocity(jj);
            }
            ++samples_written;
        }
        // Advance the dynamics
        particles.RK3(dt);
    }

    // Stop the timer
    auto end_time {std::chrono::system_clock::now()};

    // Sample at the last time
    for (int ii = 0; ii < n_particles; ++ii)
    {
        positions_samples[samples_written*n_particles +ii ] = particles.get_position(ii);
        velocity_samples[ samples_written*n_particles +ii ] = particles.get_velocity(ii);
    }
    ++samples_written;
    std::cout << "Samples Written: " << samples_written << std::endl;

    // Report computational time
    std::chrono::duration<double> computation_time {end_time - start_time};
    std::cout << "Elapsed computational time: " << computation_time.count() << std::endl;

    // Report final position and velocity
    std::cout << "The system is currently at time " << total_steps * dt << std::endl;
    std::cout << "The last particle has vel " << particles.get_velocity((n_particles -1)) << std::endl;
    std::cout << "The last particle has pos " << particles.get_position((n_particles -1)) << std::endl;

    // Save samples to binary files.
    std::ofstream sample_and_particle_size_buffer("sample_and_particle_size.bin", std::ios::out | std::ios::binary);
    std::ofstream sample_vel_file_buffer("sample_velocities.bin", std::ios::out | std::ios::binary );
    std::ofstream sample_pos_file_buffer("sample_positions.bin", std::ios::out | std::ios::binary );
    std::ofstream temp_array_buffer("initial_bath_temps.bin", std::ios::out | std::ios::binary );
    
    if(!sample_vel_file_buffer 
        || !sample_pos_file_buffer 
        || !sample_and_particle_size_buffer 
        || !temp_array_buffer )
    {
        std::cerr << "Error: Cannot create the binary files" << std::endl;
        return 1;
    }

    // Prepare array 
    int sample_and_particle_size_array[2] {n_particles, num_sample_times};
    double temp_array[2] {left_temp, right_temp};

    sample_vel_file_buffer.write( reinterpret_cast<char *>(velocity_samples), sizeof(double) * num_sample_times * n_particles );
    sample_pos_file_buffer.write( reinterpret_cast<char *>(positions_samples), sizeof(double)* num_sample_times * n_particles );
    sample_and_particle_size_buffer.write( reinterpret_cast<char *>(sample_and_particle_size_array), sizeof(int)*2);
    temp_array_buffer.write( reinterpret_cast<char *>(temp_array), sizeof(double)*2 );

    sample_vel_file_buffer.close();
    sample_pos_file_buffer.close();
    sample_and_particle_size_buffer.close();
    temp_array_buffer.close();

    delete[] positions_samples;
    delete[] velocity_samples;
    positions_samples = nullptr;
    velocity_samples = nullptr;

    return 0;
}