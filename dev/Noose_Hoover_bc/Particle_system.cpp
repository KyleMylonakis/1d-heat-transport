#include <random>
#include <iostream>
#include <cstdlib>
#include "Particle_system.h"

// Initialize and allocate the particle system
Particle_system::Particle_system(int n_particles,
    double left_temp,
    double left_NH_coupling_const,
    double right_temp,
    double right_NH_coupling_const) :
    m_n_particles {n_particles},
    m_left_temp {left_temp},
    m_left_NH_coupling_const {left_NH_coupling_const},
    m_left_NH_dynamical_var {0.0},
    m_right_temp {right_temp},
    m_right_NH_coupling_const {right_NH_coupling_const},
    m_right_NH_dynamical_var {0.0},
    m_temp_left_NH_dynamical_var {0.0},
    m_temp_right_NH_dynamical_var{0.0},
    m_rk_temp_left_NH_dynamical_var {0.0},
    m_rk_temp_right_NH_dynamical_var {0.0},
    m_k1_left_NH_dynamical_var {0.0},
    m_k1_right_NH_dynamical_var {0.0},
    m_k2_left_NH_dynamical_var {0.0},
    m_k2_right_NH_dynamical_var {0.0},
    m_k3_left_NH_dynamical_var {0.0},
    m_k3_right_NH_dynamical_var {0.0}
    {
        // Allocate the positions and velocities

        std::cout << "Allocating memory for displacements and velocities" << std::endl;
        m_positions = new double[m_n_particles];
        m_velocities = new double[m_n_particles];
        
        // Randomize the initial positions and velocities
        //std::random_device rd;
        //std::mt19937 generator(rd());
        std::mt19937 generator(10); // Fixed seed for debugging

        const double initial_avg_temp {(m_left_temp + m_right_temp)/2.0};
        const double mean {0.0};
        const double variance { initial_avg_temp };
        const double max_rand_radius { std::min(initial_avg_temp, m_a0/4.0)};

        std::normal_distribution<double> norm_dist(mean, sqrt(variance)) ;
        std::uniform_real_distribution<double> unif_dist(-max_rand_radius, max_rand_radius); 

        std::cout << "Randomizing initial displacements and velocities" << std::endl;
        for (int ii = 0; ii < m_n_particles; ++ii)
        {
            m_positions[ii] = unif_dist(generator);
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
        std::cerr << "DIVIDING BY ZerO " << std::endl;
        std::exit(1);
    }
     return 12*pow(x,-13) - 6*pow(x,-7) - 4*pow(x - m_a0,3);
    //return x - m_a0;
}

void Particle_system::f()
{
    // Store current positions, vels, and NH dynamical variables
    // in temporary arrays
    #pragma omp parallel for
        for (int ii = 0; ii < m_n_particles; ++ii)
        {
            m_temp_pos[ii] = m_positions[ii];
            m_temp_vel[ii] = m_velocities[ii];
        }
    m_temp_left_NH_dynamical_var = m_left_NH_dynamical_var;
    m_temp_right_NH_dynamical_var = m_right_NH_dynamical_var;

    // Update the new values according to f
    // only temp variables should be on the RHS
    m_velocities[0] = -Particle_system::F(m_temp_pos[1] - m_temp_pos[0] + m_a0)
        - m_temp_left_NH_dynamical_var * m_temp_vel[0];
    m_positions[0] = m_temp_vel[0];

    #pragma omp parallel for
        for (int ii = 1; ii < m_n_particles -1; ++ii)
        {
            m_velocities[ii] = Particle_system::F(m_temp_pos[ii] - m_temp_pos[ii -1] + m_a0)
                - Particle_system::F(m_temp_pos[ii + 1] - m_temp_pos[ii] + m_a0);
            m_positions[ii] = m_temp_vel[ii];
        }

    m_velocities[m_n_particles -1] = Particle_system::F(m_temp_pos[m_n_particles -1] - m_temp_pos[m_n_particles -2] + m_a0)
        - m_temp_right_NH_dynamical_var * m_temp_vel[m_n_particles -1] ;
    m_positions[m_n_particles -1] = m_temp_vel[m_n_particles -1];

    m_left_NH_dynamical_var = pow(m_left_NH_coupling_const,-2)*( (pow(m_temp_vel[0],2)/m_left_temp) - 1);
    m_right_NH_dynamical_var = pow(m_right_NH_coupling_const,-2)*( (pow(m_temp_vel[m_n_particles -1],2)/m_right_temp) - 1 );

}

// Deep copy one set of arrays to another
void Particle_system::store_pos_vel_dyn_vars(
        double* source_vel,
        double* source_pos,
        double source_left_dyn_var,
        double source_right_dyn_var,
        double* target_vel,
        double* target_pos,
        double* target_left_dyn_var,
        double* target_right_dyn_var
    )
    {
        #pragma omp parallel for
            for (int ii = 0; ii < m_n_particles; ++ii)
            {
                target_vel[ii] = source_vel[ii];
                target_pos[ii] = source_pos[ii];
            }
        *target_left_dyn_var = source_left_dyn_var;
        *target_right_dyn_var = source_right_dyn_var;
    }

void Particle_system::RK3(double h)
{
    // Store original system state
    Particle_system::store_pos_vel_dyn_vars(m_velocities, m_positions, m_left_NH_dynamical_var, m_right_NH_dynamical_var,
        m_rk_temp_vel, m_rk_temp_pos, &m_rk_temp_left_NH_dynamical_var, &m_rk_temp_right_NH_dynamical_var);

    
    // k1 = f(y)
    Particle_system::f();
    Particle_system::store_pos_vel_dyn_vars(m_velocities, m_positions, m_left_NH_dynamical_var, m_right_NH_dynamical_var,
        m_k1_vel, m_k1_pos, &m_k1_left_NH_dynamical_var, &m_k1_right_NH_dynamical_var);

    // update current state to y + h/2 * k1
    #pragma omp parallel for
        for (int ii = 0; ii < m_n_particles; ++ii)
        {
            m_velocities[ii] = m_rk_temp_vel[ii] + (h/2.0) * m_k1_vel[ii];
            m_positions[ii] = m_rk_temp_pos[ii] + (h/2.0) * m_k1_pos[ii];
        }
    m_left_NH_dynamical_var = m_rk_temp_left_NH_dynamical_var + (h/2.0) * m_k1_left_NH_dynamical_var;
    m_right_NH_dynamical_var = m_rk_temp_right_NH_dynamical_var + (h/2.0) * m_k1_right_NH_dynamical_var;

    // Compute f(y + h/2 * k1) and store in k2
    Particle_system::f();
    Particle_system::store_pos_vel_dyn_vars(m_velocities, m_positions, m_left_NH_dynamical_var, m_right_NH_dynamical_var,
        m_k2_vel, m_k2_pos, &m_k2_left_NH_dynamical_var, &m_k2_right_NH_dynamical_var);

    // Update current state to y + h(-k1 + 2k2)
    #pragma omp parallel for
        for (int ii = 0; ii < m_n_particles; ++ii)
        {
            m_velocities[ii] = m_rk_temp_vel[ii] + h * (- m_k1_vel[ii] + 2 * m_k2_vel[ii]);
            m_positions[ii] = m_rk_temp_pos[ii] + h * (- m_k1_pos[ii] + 2 * m_k2_pos[ii]);
        }
    m_left_NH_dynamical_var = m_rk_temp_left_NH_dynamical_var + h * (-m_k1_left_NH_dynamical_var + 2 * m_k2_left_NH_dynamical_var);
    m_right_NH_dynamical_var = m_rk_temp_right_NH_dynamical_var + h * (-m_k1_right_NH_dynamical_var + 2 * m_k2_right_NH_dynamical_var);

    // Compute f(y + h(-k1 + 2k2) ) and store in k3
    Particle_system::f();
    Particle_system::store_pos_vel_dyn_vars(m_velocities, m_positions, m_left_NH_dynamical_var, m_right_NH_dynamical_var,
        m_k3_vel, m_k3_pos, &m_k3_left_NH_dynamical_var, &m_k3_right_NH_dynamical_var);
    
    // Compute the final update at the next time step
    #pragma omp parallel for
        for (int ii = 0; ii < m_n_particles; ++ii)
        {
            m_velocities[ii] = m_rk_temp_vel[ii]
                + h * ( (1.0/6.0)*m_k1_vel[ii] + (2.0/3.0)*m_k2_vel[ii] + (1.0/6.0)*m_k3_vel[ii] );
            m_positions[ii] = m_rk_temp_pos[ii]
                + h * ( (1.0/6.0)*m_k1_pos[ii] + (2.0/3.0)*m_k2_pos[ii] + (1.0/6.0)*m_k3_pos[ii] );
        }
    m_left_NH_dynamical_var = m_rk_temp_left_NH_dynamical_var
        + h * ( (1.0/6.0)*m_k1_left_NH_dynamical_var + (2.0/3.0)*m_k2_left_NH_dynamical_var + (1.0/6.0)*m_k3_left_NH_dynamical_var );
    m_right_NH_dynamical_var = m_rk_temp_right_NH_dynamical_var
        + h * ( (1.0/6.0)*m_k1_right_NH_dynamical_var + (2.0/3.0)*m_k2_right_NH_dynamical_var + (1.0/6.0)*m_k3_right_NH_dynamical_var );

}

double Particle_system::get_position(int index) const
{
    return m_positions[index];
}

double Particle_system::get_velocity(int index) const
{
    return m_velocities[index];
}

double Particle_system::compute_mean(double* array, int length)
{
    double mean {0.0};
    for (int ii = 0; ii < length; ++ii)
    {
        mean += array[ii];
    }
    return mean/length;
}

double Particle_system::compute_variance(double* array, int length)
{
    double var {0.0};
    double mean {compute_mean(array, length)};
    for (int ii = 0; ii < length; ++ii)
    {
        var += std::pow( array[ii] - mean,2)    ;
    }
    return var/length;
}