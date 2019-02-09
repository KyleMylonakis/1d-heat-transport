/*
    * The class defining the particle system
    * This class is responsible for initializing
    * the system, defining the potential, defining
    * the interaction between particles, and how the system
    * advances the dynamics
*/


#include <iostream>
#include <cstdlib>
#include "Particle_system.h"

// Initialize and allocate the particle system
Particle_system::Particle_system(int n_particles,
    double left_temp,
    double right_temp,
    double langevin_coupling_const,
    std::string boundary_condition) :
    m_n_particles {n_particles},
    m_left_temp {left_temp},
    m_right_temp {right_temp},
    m_langevin_coupling_const { langevin_coupling_const },
    m_left_langevin_variance { 2 * m_langevin_coupling_const * m_left_temp },
    m_right_langevin_variance  {2 * m_langevin_coupling_const * m_right_temp }
    {
        // Allocate the positions and velocities

        std::cout << "Allocating memory for displacements and velocities" << std::endl;
        m_positions = new double[m_n_particles]();
        m_velocities = new double[m_n_particles]();
        
        // Randomize the initial positions and velocities
        std::random_device rd;
        std::mt19937 generator(rd());
        //std::mt19937 generator(10); // Fixed seed for debugging

        const double initial_avg_temp {(m_left_temp + m_right_temp)/2.0};
        const double mean {0.0};
        const double variance { initial_avg_temp };
        const double max_rand_radius { std::min(initial_avg_temp, m_a0/4.0)};
        

        std::normal_distribution<double> norm_dist(mean, sqrt(variance)) ;
        std::uniform_real_distribution<double> unif_dist(-max_rand_radius, max_rand_radius); 

        std::cout << "Randomizing initial displacements and velocities" << std::endl;
        for (int ii = 0; ii < m_n_particles; ++ii)
        {
            //m_positions[ii] = unif_dist(generator);
            m_velocities[ii] = norm_dist(generator);
        }

        // Allocate memory and initialize vars for f() method
        std::cout << "Allocating memory for f() method" << std::endl;
        m_temp_vel = new double[m_n_particles]();
        m_temp_pos = new double[m_n_particles]();

        // Allocate memory and initialize vars for RK3 method
        m_rk_temp_vel = new double[m_n_particles]();
        m_rk_temp_pos = new double[m_n_particles]();
        m_k1_vel = new double[m_n_particles]();
        m_k1_pos = new double[m_n_particles]();
        m_k2_vel = new double[m_n_particles]();
        m_k2_pos = new double[m_n_particles]();
        m_k3_vel = new double[m_n_particles]();
        m_k3_pos = new double[m_n_particles]();

        // Set internal_variables for the boundary condition
        if (boundary_condition == "fixed")
        {
            m_boundary_condition = 1;
        }
        else if (boundary_condition == "free")
        {
            m_boundary_condition = 0;
        }
        else
        {
            std::cerr << "ERROR: Invalid choice of boundary condition. Exiting" << std::endl;
            Particle_system::~Particle_system();
            std::exit(1);
        }
        
        m_left_stochastic_bath = std::normal_distribution<double>(mean, sqrt(m_left_langevin_variance));
        m_right_stochastic_bath = std::normal_distribution<double>(mean, sqrt(m_right_langevin_variance));
        m_generator = std::mt19937(m_rd());
        //m_generator = std::mt19937(10);
    }

Particle_system::~Particle_system()
{
    // Deallocate the positions and velocities
    std::cout << "Deallocating positions and velocities" << std::endl;
    delete[] m_positions;
    m_positions = nullptr;
    delete[] m_velocities;
    m_velocities = nullptr;

    // Deallocate the temporary variables for the f() method
    std::cout << "Deallocate the temporary variables for the f() method" << std::endl;
    delete[] m_temp_pos;
    m_temp_pos = nullptr;
    delete[] m_temp_vel;
    m_temp_vel = nullptr;

    // Deallocate the temporary variable for the RK3 method
    delete[] m_rk_temp_vel;
    delete[] m_rk_temp_pos;
    delete[] m_k1_vel;
    delete[] m_k1_pos;
    delete[] m_k2_vel;
    delete[] m_k2_pos;
    delete[] m_k3_vel;
    delete[] m_k3_pos;

    m_rk_temp_vel = nullptr;
    m_rk_temp_pos = nullptr;
    m_k1_vel = nullptr;
    m_k1_pos = nullptr;
    m_k2_vel = nullptr;
    m_k2_pos = nullptr;
    m_k3_vel = nullptr;
    m_k3_pos = nullptr;
}

// Interparticle Force: Comes from LJ potential with confinement
// F = - grad(U)
double Particle_system::F(double x) const
{
    
    if (!std::isfinite(12*pow(x,-13) - 6*pow(x,-7) - 4*pow(x - m_a0,3) ))
    {
        std::cerr << "DIVIDING BY ZERO " << std::endl;
        std::exit(1);
    }
    return 12*pow(x,-13) - 6*pow(x,-7) - 4*pow(x - m_a0,3);
    
   /*
    double out;
    out = 6*std::pow(x,-7);
    out *= (2*std::pow(x,-6) - 1 );
    out += 4*std::pow(x - m_a0,3);
    if (!std::isfinite(out))
    {
        std::cerr << "DIVIDING BY ZERO " << std::endl;
        std::exit(1);
    }
    return out;
    */
    //return -x + m_a0;
}

void Particle_system::f(double h)
{
    // Store current positions, vels, and NH dynamical variables
    // in temporary arrays
    #pragma omp parallel for
        for (int ii = 0; ii < m_n_particles; ++ii)
        {
            m_temp_pos[ii] = m_positions[ii];
            m_temp_vel[ii] = m_velocities[ii];
        }
    
    // Update the new values according to f
    // only temp variables should be on the RHS
    m_velocities[0] = m_boundary_condition * Particle_system::F(m_temp_pos[0] + m_a0) 
        - Particle_system::F(m_temp_pos[1] - m_temp_pos[0] + m_a0)
        //- m_langevin_coupling_const * m_temp_vel[0] // Commented out for IMEX
        + m_left_stochastic_bath(m_generator)/sqrt(h);
    m_positions[0] = m_temp_vel[0];

    #pragma omp parallel for
        for (int ii = 1; ii < m_n_particles -1; ++ii)
        {
            m_velocities[ii] = Particle_system::F(m_temp_pos[ii] - m_temp_pos[ii -1] + m_a0)
                - Particle_system::F(m_temp_pos[ii + 1] - m_temp_pos[ii] + m_a0);
            m_positions[ii] = m_temp_vel[ii];
        }

    m_velocities[m_n_particles -1] = Particle_system::F(m_temp_pos[m_n_particles -1] - m_temp_pos[m_n_particles -2] + m_a0)
        - m_boundary_condition * Particle_system::F( - m_temp_pos[m_n_particles - 1] + m_a0)
        //- m_langevin_coupling_const * m_temp_vel[m_n_particles -1] // Commented out for IMEX
        + m_right_stochastic_bath(m_generator)/sqrt(h);
    m_positions[m_n_particles -1] = m_temp_vel[m_n_particles -1];

}

// Deep copy one set of arrays to another
void Particle_system::store_pos_vel_dyn_vars(
        double* source_vel,
        double* source_pos,
        double* target_vel,
        double* target_pos
    )
    {
        #pragma omp parallel for
            for (int ii = 0; ii < m_n_particles; ++ii)
            {
                target_vel[ii] = source_vel[ii];
                target_pos[ii] = source_pos[ii];
            }
    }

void Particle_system::RK3(double h)
{
    std::cerr << "DO NOT USE WITHOUT FIXING f() " << std::endl;
    exit(1);
    // Store original system state
    Particle_system::store_pos_vel_dyn_vars(m_velocities, m_positions,
        m_rk_temp_vel, m_rk_temp_pos);

    
    // k1 = f(y)
    Particle_system::f(h);
    Particle_system::store_pos_vel_dyn_vars(m_velocities, m_positions,
        m_k1_vel, m_k1_pos);

    // update current state to y + h/2 * k1
    #pragma omp parallel for
        for (int ii = 0; ii < m_n_particles; ++ii)
        {
            m_velocities[ii] = m_rk_temp_vel[ii] + (h/2.0) * m_k1_vel[ii];
            m_positions[ii] = m_rk_temp_pos[ii] + (h/2.0) * m_k1_pos[ii];
        }

    // Compute f(y + h/2 * k1) and store in k2
    Particle_system::f(h);
    Particle_system::store_pos_vel_dyn_vars(m_velocities, m_positions,
        m_k2_vel, m_k2_pos);

    // Update current state to y + h(-k1 + 2k2)
    #pragma omp parallel for
        for (int ii = 0; ii < m_n_particles; ++ii)
        {
            m_velocities[ii] = m_rk_temp_vel[ii] + h * (- m_k1_vel[ii] + 2 * m_k2_vel[ii]);
            m_positions[ii] = m_rk_temp_pos[ii] + h * (- m_k1_pos[ii] + 2 * m_k2_pos[ii]);
        }

    // Compute f(y + h(-k1 + 2k2) ) and store in k3
    Particle_system::f(h);
    Particle_system::store_pos_vel_dyn_vars(m_velocities, m_positions,
        m_k3_vel, m_k3_pos);
    
    // Compute the final update at the next time step
    #pragma omp parallel for
        for (int ii = 0; ii < m_n_particles; ++ii)
        {
            m_velocities[ii] = m_rk_temp_vel[ii]
                + h * ( (1.0/6.0)*m_k1_vel[ii] + (2.0/3.0)*m_k2_vel[ii] + (1.0/6.0)*m_k3_vel[ii] );
            m_positions[ii] = m_rk_temp_pos[ii]
                + h * ( (1.0/6.0)*m_k1_pos[ii] + (2.0/3.0)*m_k2_pos[ii] + (1.0/6.0)*m_k3_pos[ii] );
        }

}

void Particle_system::milstein(double h)
{
    // Store the current position and velocity
    Particle_system::store_pos_vel_dyn_vars(m_velocities, m_positions,
        m_rk_temp_vel, m_rk_temp_pos);

    // Computes (a*dt + b*dW)/h
    Particle_system::f(h);

    // Computes the update based on the Milstein scheme.
    // First and last velocities are using IMEX scheme
    m_velocities[0] = (m_rk_temp_vel[0] + h * m_velocities[0])/(1.0 + m_langevin_coupling_const*h);
    m_positions[0] = m_rk_temp_pos[0] + h * m_positions[0];
    #pragma omp parallel for
        for (int ii = 1; ii < m_n_particles - 1; ++ii)
        {
            m_velocities[ii] = m_rk_temp_vel[ii] + h * m_velocities[ii];
            m_positions[ii] = m_rk_temp_pos[ii] + h * m_positions[ii];
        }
    m_velocities[m_n_particles - 1] = (m_rk_temp_vel[m_n_particles - 1] + h * m_velocities[m_n_particles - 1])/(1.0 + m_langevin_coupling_const*h);
    m_positions[m_n_particles -1] = m_rk_temp_pos[m_n_particles -1] + h * m_positions[m_n_particles -1];
}

double Particle_system::get_position(int index) const
{
    return m_positions[index];
}

double Particle_system::get_velocity(int index) const
{
    return m_velocities[index];
}

double Particle_system::compute_mean(double* array, int length) const
{
    double mean {0.0};
    for (int ii = 0; ii < length; ++ii)
    {
        mean += array[ii];
    }
    return mean/length;
}

double Particle_system::compute_variance(double* array, int length) const
{
    double var {0.0};
    double mean {compute_mean(array, length)};
    for (int ii = 0; ii < length; ++ii)
    {
        var += std::pow( array[ii] - mean,2)    ;
    }
    return var/length;
}

double Particle_system::sample_left_bath()
{
    return m_left_stochastic_bath(m_generator);
}

double Particle_system::sample_right_bath()
{
    return m_right_stochastic_bath(m_generator);
}