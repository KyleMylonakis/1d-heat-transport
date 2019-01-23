#pragma once

#include <fstream>
#include <iostream>
#include <string>
#include <cmath>
#include <limits>
#include <iomanip>
#include <string>
#include "Weights_and_Nodes.h"


// Initialize the constants
extern const double g_dt {0.001}; // Time step of the simulation
const int dt_inv {static_cast<int>(1/g_dt)};      // Must be 1/dt as an int
const int t_max {100};        // Final time of the simulation
extern const double g_k_sqrd {156*pow(2,-7.0/3.0) - 42*pow(2, -4.0/3.0) }; // Coming from the second derivative of the Lennard Jones 
extern const double g_k {sqrt(g_k_sqrd)}; // Coming from the second derivative of the Lennard Jones 
extern const double a0 {pow(2,1.0/6.0)}; // Equilibrium spacing
//const int n_particles {100};      // The number of particles in the resolved system // YOU NEED TO CHANGE THIS OVER IN FUNCTIONS.CPP AS WELL
//extern const int n_particles {100};      // The number of particles in the resolved system // YOU NEED TO CHANGE THIS OVER IN FUNCTIONS.CPP AS WELL
extern const int n_particles {1024};      // FOR DEBUGGING ONLY
extern const long long int total_steps {dt_inv*t_max + 1}; // Number of steps in the simulation
int const t_equi {t_max/2};          // Time until "equilibriation". The time when I start sampling the data
long long int const steps_to_equi {t_equi *dt_inv + 1}; // Number of steps until "equilibriation"
int const sample_freq {dt_inv/10};      // Sample frequency after equilibriation
long long int const total_samples { (t_max -t_equi)*sample_freq }; // Total number of samples
int const n_int_particles {n_particles -2};
//const int n_bath_particles {9000}; // Number of particles being sampled in the baths. Must be the same as in functions.cpp as well
extern const long long int n_bath_particles {90000}; // Number of particles being sampled in the baths. Must be the same as in functions.cpp as well
//extern const long long int n_bath_particles {100}; // FOR DEBUGGING ONLY
extern const long long int max_trap_steps {total_steps}; // Max number of steps to take for the trapezoidal rule
double trunc_beta_norm {0}; // The value of the integral of beta from 0 to dt*max_trap_steps
double trunc_theta_norm {0}; // The value of the integral of theta from 0 to dt*max_trap_steps
int quarter_grid_steps {9};
int half_grid_steps {quarter_grid_steps/2 + 1}; // This uses int division


// Allocate the Arrays
double* l_bath_init_pos = new double[n_bath_particles]; // Initial displacements of particles int the left bath
double* l_bath_init_vel = new double[n_bath_particles]; // Initial velocities of particles int the left bath
double* r_bath_init_pos = new double[n_bath_particles]; // Initial displacements of particles int the right bath
double* r_bath_init_vel = new double[n_bath_particles]; // Initial velocities of particles int the right bath

double* first_vel = new double[total_steps];
double* last_vel = new double[total_steps];
double* positions = new double[n_particles];
double* int_velocities = new double[n_int_particles];
double* vel_samples = new double[total_samples*n_particles];
double* pos_samples = new double[total_samples*n_particles];

double* forces = new double[n_particles];
double* beta = new double[total_steps];
double* shifted_theta = new double[total_steps];
double* theta = new double[total_steps]();

double* a_0 = new double[2*total_steps-1]();
double* b_0 = new double[2*total_steps-1]();
double* a_bath = new double[n_bath_particles -1]();
double* b_bath = new double[n_bath_particles -1]();


double* quarter_grid_first_vel = new double[quarter_grid_steps]();
double* quarter_grid_last_vel = new double[quarter_grid_steps]();
double* quarter_grid_int_vel = new double[n_int_particles]();
double* quarter_grid_pos = new double[n_particles]();
double* store_quarter_grid_int_vel = new double[n_int_particles * quarter_grid_steps]();
double* store_quarter_grid_pos = new double[n_particles * quarter_grid_steps]();

double* half_grid_first_vel = new double[half_grid_steps]();
double* half_grid_last_vel = new double[half_grid_steps]();
double* half_grid_int_vel = new double[n_int_particles]();
double* half_grid_pos = new double[n_particles]();
double* store_half_grid_int_vel = new double[n_int_particles * half_grid_steps]();
double* store_half_grid_pos = new double[n_particles * half_grid_steps]();




#include "function_forward_declarations.h"


// Function to cleanup all the dynamic arrays
void cleanup(void)
{
    delete[] l_bath_init_pos;
    delete[] l_bath_init_vel;
    delete[] r_bath_init_pos;
    delete[] r_bath_init_vel;
    delete[] first_vel;
    delete[] last_vel;
    delete[] int_velocities;
    delete[] positions;
    delete[] vel_samples;
    delete[] pos_samples;
    delete[] forces;
    delete[] beta;
    delete[] theta;
    delete[] shifted_theta;
    delete[] a_0;
    delete[] b_0;
    delete[] a_bath;
    delete[] b_bath;
    delete[] quarter_grid_first_vel;
    delete[] quarter_grid_last_vel;
    delete[] quarter_grid_int_vel;
    delete[] quarter_grid_pos;
    delete[] half_grid_first_vel;
    delete[] half_grid_last_vel;
    delete[] half_grid_int_vel;
    delete[] half_grid_pos;
    delete[] store_quarter_grid_int_vel;
    delete[] store_quarter_grid_pos;
    delete[] store_half_grid_int_vel;
    delete[] store_half_grid_pos;
    

    l_bath_init_pos = nullptr;
    l_bath_init_vel = nullptr;
    r_bath_init_pos = nullptr;
    r_bath_init_vel = nullptr;
    first_vel = nullptr;
    last_vel = nullptr;
    int_velocities = nullptr;
    positions = nullptr;
    vel_samples = nullptr;
    pos_samples = nullptr;
    forces = nullptr;
    beta = nullptr;
    theta = nullptr;
    shifted_theta = nullptr;
    a_0 = nullptr;
    b_0 = nullptr;
    a_bath = nullptr;
    b_bath = nullptr;
    quarter_grid_first_vel = nullptr;
    quarter_grid_last_vel = nullptr;
    quarter_grid_int_vel = nullptr;
    quarter_grid_pos = nullptr;
    half_grid_first_vel = nullptr;
    half_grid_last_vel = nullptr;
    half_grid_int_vel = nullptr;
    half_grid_pos = nullptr;
    store_quarter_grid_int_vel = nullptr;
    store_quarter_grid_pos = nullptr;
    store_half_grid_int_vel = nullptr;
    store_half_grid_pos = nullptr;
}
