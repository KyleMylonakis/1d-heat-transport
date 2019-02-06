#pragma once

#include <cmath>

class Particle_system
{
private:
    int m_n_particles;
    double* m_positions;
    double* m_velocities;

    double m_left_NH_dynamical_var;
    double m_left_NH_coupling_const;
    double m_left_temp;

    double m_right_NH_dynamical_var;
    double m_right_NH_coupling_const;
    double m_right_temp;

    const double m_a0 {pow(2,1.0/6.0)};

    // Extra arrays for f() method
    double* m_temp_vel;
    double* m_temp_pos;
    double m_temp_left_NH_dynamical_var;
    double m_temp_right_NH_dynamical_var;

    // Temporary arrays & vars for RK3 method
    double* m_rk_temp_vel;
    double* m_rk_temp_pos;
    double m_rk_temp_left_NH_dynamical_var;
    double m_rk_temp_right_NH_dynamical_var;

    double* m_k1_vel;
    double* m_k1_pos;
    double m_k1_left_NH_dynamical_var;
    double m_k1_right_NH_dynamical_var;

    double* m_k2_vel;
    double* m_k2_pos;
    double m_k2_left_NH_dynamical_var;
    double m_k2_right_NH_dynamical_var;

    double* m_k3_vel;
    double* m_k3_pos;
    double m_k3_left_NH_dynamical_var;
    double m_k3_right_NH_dynamical_var;

    // Interparticle Force: Comes from LJ potential with confinement
    // F = - grad(U)
    double F(double x) const;

    // Our ODE has the vector form y' = f(y)
    // f(y) will be a function which changes the state of the object
    // from y to f(y).
    void f();


    // Deep copy one set of arrays to another
    void store_pos_vel_dyn_vars(
        double* source_vel,
        double* source_pos,
        double source_left_dyn_var,
        double source_right_dyn_var,
        double* target_vel,
        double* target_pos,
        double* target_left_dyn_var,
        double* target_right_dyn_var
    );

public:

    // Delete the default constructor to require initialization 
    // with all input parameters
    Particle_system() = delete;

    // Initialize and allocate the particle system
    Particle_system(int n_particles,
        double left_temp,
        double left_NH_coupling_const,
        double right_temp,
        double right_NH_coupling_const);

    ~Particle_system();

    void RK3(double h);

    double get_position(int index) const;
    double get_velocity(int index) const;

    double compute_mean(double* array, int length);
    double compute_variance(double* array, int length);

};