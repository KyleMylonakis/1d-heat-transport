#pragma once

#include <cmath>
#include <string>
#include <random>

class Particle_system
{
private:
    int m_n_particles;
    double* m_positions;
    double* m_velocities;

    const double m_left_temp;
    const double m_right_temp;
    const double m_langevin_coupling_const;

    const double m_a0 {pow(2,1.0/6.0)};
    const double m_nu {10.0};
    const double m_c0;
    const double m_c1;
    const double m_c2;
    
    // Langevin parameters and noise
    const double m_left_langevin_variance ;
    const double m_right_langevin_variance;
    std::normal_distribution<double> m_left_stochastic_bath;
    std::normal_distribution<double> m_right_stochastic_bath;
    std::random_device m_rd;
    std::mt19937 m_generator;

    double m_time_step;
    
    int m_damping;

    // Extra arrays for f() method
    double* m_temp_vel;
    double* m_temp_pos;
    
    // Temporary arrays & vars for RK3 method
    double* m_rk_temp_vel;
    double* m_rk_temp_pos;   
    double* m_k1_vel;
    double* m_k1_pos;
    double* m_k2_vel;
    double* m_k2_pos;
    double* m_k3_vel;
    double* m_k3_pos;

    // Arrays for Velocity Verlet with Damping
    double* m_accel;
    double* m_forces;
    double* m_forces_prev_step;
    
    // Integer representing the fixed or free boundary condition
    int m_boundary_condition;

    // Interparticle Force: Comes from LJ potential with confinement
    // F = - grad(U)
    double F(double x) const;

    // Our ODE has the vector form y' = f(y)
    // f(y) will be a function which changes the state of the object
    // from y to f(y).
    void f(double h);

    // accel = forces - damping
    void compute_accel(double h);
    void compute_forces(double h);
    
    // Deep copy
    void deep_copy(double* source, double* target, int length);

    // Deep copy one set of arrays to another
    void store_pos_vel_dyn_vars(
        double* source_vel,
        double* source_pos,
        double* target_vel,
        double* target_pos
    );

public:

    // Delete the default constructor to require initialization 
    // with all input parameters
    Particle_system() = delete;

    // Initialize and allocate the particle system
    Particle_system(int n_particles,
        double left_temp,
        double right_temp,
        double langevin_coupling_const,
        std::string boundary_condition,
        double time_step,
        std::string damping);

    ~Particle_system();

    void RK3(double h);
    void milstein (double h);
    void velocity_verlet_with_damping();

    double get_position(int index) const;
    double get_velocity(int index) const;

    double compute_mean(double* array, int length) const;
    double compute_variance(double* array, int length) const;
    double sample_left_bath();
    double sample_right_bath();

    void velocity_verlet_integrating_factor();

};