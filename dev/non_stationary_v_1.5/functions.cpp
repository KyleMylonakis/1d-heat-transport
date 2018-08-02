/*
* This file contains all the auxillary functions for the RK3 integrators
*/

//#include "constant_forward_declarations.h"

#include <iostream>
#include <cmath>
#include <algorithm>
#include <limits>
#include <omp.h>
#include <iomanip>

//#define NUM_THREADS 1

extern const double g_dt; // Time step
extern const double g_k_sqrd; // Coming from the second derivative of the Lennard Jones 
extern const double g_k;
extern const double a0; // Equilibrium spacing
extern const int n_particles;      // The number of particles in the resolved system
const int n_int_particles = n_particles - 2;
extern const int n_bath_particles;
extern const long long int max_trap_steps;
extern const long long int total_steps;



// Define the function theta
double theta_f(double t)
{
    return (std::cyl_bessel_j(1,2*g_k*t))/(g_k*t);
}

double normalize_theta(double *theta, double scale_dt)
{
    if (max_trap_steps >= total_steps )
    {
        return 1.0;
    }
    else
    {
        double output { 0 };
        for ( long long int j = 1; j < std::min(max_trap_steps,total_steps) ; ++j ) 
            {
                output += theta[j];
            }
        output += ( 0.5 * theta[0] ) + ( 0.5 * theta[max_trap_steps] );
        output = scale_dt*g_dt*output;
        output = output*g_k;
        return output;
    }
    
}

// Define the function beta
double beta_f(double t)
{
    //if (t > 0 )
    //{
        return 2*( std::cyl_bessel_j(1,2*g_k*t)/(g_k*t*t) - std::cyl_bessel_j(0,2*g_k*t)/t);   
    //}
    //else 
    //{
    //    return 0.0;
    //}
}


// Define the function a0(t)
double a_0_f(double t)
{
    return beta_f(t)/g_k_sqrd;
}

// Define the function b0(t) (which is the derivative of a0(t))
double b_0_f(double t)
{
    return ( 4*std::cyl_bessel_j(1,2*g_k*t)/(g_k*t) 
            + 6*std::cyl_bessel_j(0,2*g_k*t)/(g_k_sqrd*t*t) 
            - 6*std::cyl_bessel_j(1,2*g_k*t)/(g_k*g_k_sqrd*t*t*t)
           );
}


// Define the vector f
void f(double* a_in, double* b_in, double t)
{
    double* temp_a = new double[n_bath_particles -1];
    double* temp_b = new double[n_bath_particles -1];

    for (long long int j = 0; j < n_bath_particles -1; ++j)
    {
        temp_a[j] = b_in[j];
    }

    temp_b[0] = g_k_sqrd*(a_in[1] - 2*a_in[0] + a_0_f(t) );
    //std::cout << temp_b[0] << std::endl;
    for (long long int j = 1; j < n_bath_particles -2; ++j)
    {
        temp_b[j] = g_k_sqrd*( a_in[j-1] -2*a_in[j] + a_in[j+1]);
    }
    temp_b[n_bath_particles -2] = -g_k*b_in[n_bath_particles -2] - g_k_sqrd*a_in[n_bath_particles -2] + g_k_sqrd*a_in[n_bath_particles -3];

    for (long long int j = 0; j < n_bath_particles -1; ++j)
    {
        a_in[j] = temp_a[j];
        b_in[j] = temp_b[j];
        //std::cout << temp_b[0] << std::endl;
        //std::cout << b_in[0] << std::endl;
    }

    delete[] temp_a;
    delete[] temp_b;
    temp_a = nullptr;
    temp_b = nullptr;
}


void f_time_zero(double* a_in, double* b_in, double t)
{
    double* temp_a = new double[n_bath_particles -1];
    double* temp_b = new double[n_bath_particles -1];

    for (long long int j = 0; j < n_bath_particles -1; ++j)
    {
        temp_a[j] = b_in[j];
    }

    temp_b[0] = g_k_sqrd*(a_in[1] - 2*a_in[0] );
    //std::cout << temp_b[0] << std::endl;
    for (long long int j = 1; j < n_bath_particles -2; ++j)
    {
        temp_b[j] = g_k_sqrd*( a_in[j-1] -2*a_in[j] + a_in[j+1]);
    }
    temp_b[n_bath_particles -2] = -g_k*b_in[n_bath_particles -2] - g_k_sqrd*a_in[n_bath_particles -2] + g_k_sqrd*a_in[n_bath_particles -3];

    for (long long int j = 0; j < n_bath_particles -1; ++j)
    {
        a_in[j] = temp_a[j];
        b_in[j] = temp_b[j];
        //std::cout << temp_b[0] << std::endl;
        //std::cout << b_in[0] << std::endl;
    }

    delete[] temp_a;
    delete[] temp_b;
    temp_a = nullptr;
    temp_b = nullptr;
}

// Perform RK3 integration for half a step
void update_bath(double* a_0, double* b_0, double* a_bath, double* b_bath, double current_time)
{
    double* temp_a = new double[n_bath_particles -1];
    double* temp_b = new double[n_bath_particles -1];
    double* alpha_a = new double[n_bath_particles -1];
    double* alpha_b = new double[n_bath_particles -1];
    double* beta_a = new double[n_bath_particles -1];
    double* beta_b = new double[n_bath_particles -1];
    
    #pragma omp parallel for
        for (long long int j = 0; j < n_bath_particles -1; ++j) {
            temp_a[j] = a_bath[j];
            temp_b[j] = b_bath[j];
        }

    // Makes a_bath and b_bath the vector k_1
    f(a_bath, b_bath, current_time);

    #pragma omp parallel for
        for (long long int j = 0; j < n_bath_particles -1; ++j) {
            alpha_a[j] = temp_a[j] + 0.5 * g_dt * 0.5 * a_bath[j];
            alpha_b[j] = temp_b[j] + 0.5 * g_dt * 0.5 * b_bath[j];
        }

    // makes alpha_a and alpha_b the vector k_2
    f(alpha_a, alpha_b, current_time + g_dt /4.0);

    #pragma omp parallel for
        for (long long int j = 0; j < n_bath_particles -1; ++j) {
            beta_a[j] = temp_a[j] - (g_dt * 0.5 * a_bath[j]) + (2 * g_dt * 0.5 * alpha_a[j]);
            beta_b[j] = temp_b[j] - (g_dt * 0.5 * b_bath[j]) + (2 * g_dt * 0.5 * alpha_b[j]);
        }

    // makes beta_a and beta_b the vector k_3
    f(beta_a,beta_b,current_time + g_dt * 0.5);

    #pragma omp parallel for
        for (long long int j = 0; j < n_bath_particles -1; ++j) {
            a_bath[j] = temp_a[j] + g_dt * 0.5 * ( (1.0/6.0)*a_bath[j] + (2.0/3.0)*alpha_a[j] + (1.0/6.0)*beta_a[j] );
            b_bath[j] = temp_b[j] + g_dt * 0.5 * ( (1.0/6.0)*b_bath[j] + (2.0/3.0)*alpha_b[j] + (1.0/6.0)*beta_b[j] );
        }



    delete[] temp_a;
    delete[] temp_b;
    delete[] alpha_a;
    delete[] alpha_b;
    delete[] beta_a;
    delete[] beta_b;
    temp_a = nullptr;
    temp_b = nullptr;
    alpha_a = nullptr;
    alpha_b = nullptr;
    beta_a = nullptr;
    beta_b = nullptr;
} 

void update_bath_time_zero(double* a_0, double* b_0, double* a_bath, double* b_bath, double current_time)
{
    double* temp_a = new double[n_bath_particles -1];
    double* temp_b = new double[n_bath_particles -1];
    double* alpha_a = new double[n_bath_particles -1];
    double* alpha_b = new double[n_bath_particles -1];
    double* beta_a = new double[n_bath_particles -1];
    double* beta_b = new double[n_bath_particles -1];
    

    for (long long int j = 0; j < n_bath_particles -1; ++j)
    {
        temp_a[j] = a_bath[j];
        temp_b[j] = b_bath[j];
    }

    // Makes a_bath and b_bath the vector k_1
    f_time_zero(a_bath, b_bath, current_time);

    for (long long int j = 0; j < n_bath_particles -1; ++j)
    {
        alpha_a[j] = temp_a[j] + 0.5 * g_dt * 0.5 * a_bath[j];
        alpha_b[j] = temp_b[j] + 0.5 * g_dt * 0.5 * b_bath[j];
    }

    // makes alpha_a and alpha_b the vector k_2
    f(alpha_a, alpha_b, current_time + (0.5 * g_dt * 0.5));

    for (long long int j = 0; j < n_bath_particles -1; ++j)
    {
        beta_a[j] = temp_a[j] - (g_dt * 0.5 * a_bath[j]) + (2 * g_dt * 0.5 * alpha_a[j]);
        beta_b[j] = temp_b[j] - (g_dt * 0.5 * b_bath[j]) + (2 * g_dt * 0.5 * alpha_b[j]);
    }

    // makes beta_a and beta_b the vector k_3
    f(beta_a,beta_b,current_time + g_dt * 0.5);

    for (long long int j = 0; j < n_bath_particles -1; ++j)
    {
        a_bath[j] = temp_a[j] + g_dt * 0.5 * ( (1.0/6.0)*a_bath[j] + (2.0/3.0)*alpha_a[j] + (1.0/6.0)*beta_a[j] );
        b_bath[j] = temp_b[j] + g_dt * 0.5 * ( (1.0/6.0)*b_bath[j] + (2.0/3.0)*alpha_b[j] + (1.0/6.0)*beta_b[j] );
    }



    delete[] temp_a;
    delete[] temp_b;
    delete[] alpha_a;
    delete[] alpha_b;
    delete[] beta_a;
    delete[] beta_b;
    temp_a = nullptr;
    temp_b = nullptr;
    alpha_a = nullptr;
    alpha_b = nullptr;
    beta_a = nullptr;
    beta_b = nullptr;
} 

// Perform RK3 integration for half a step
void update_bath_quarter(double* a_0, double* b_0, double* a_bath, double* b_bath, double current_time)
{
    double* temp_a = new double[n_bath_particles -1];
    double* temp_b = new double[n_bath_particles -1];
    double* alpha_a = new double[n_bath_particles -1];
    double* alpha_b = new double[n_bath_particles -1];
    double* beta_a = new double[n_bath_particles -1];
    double* beta_b = new double[n_bath_particles -1];
    

    for (long long int j = 0; j < n_bath_particles -1; ++j)
    {
        temp_a[j] = a_bath[j];
        temp_b[j] = b_bath[j];
    }

    // Makes a_bath and b_bath the vector k_1
    f(a_bath, b_bath, current_time);

    for (long long int j = 0; j < n_bath_particles -1; ++j)
    {
        alpha_a[j] = temp_a[j] + (0.5 * g_dt * 0.25 * a_bath[j]);
        alpha_b[j] = temp_b[j] + (0.5 * g_dt * 0.25 * b_bath[j]);
    }

    // makes alpha_a and alpha_b the vector k_2
    f(alpha_a, alpha_b, current_time + (g_dt /8.0));

    for (long long int j = 0; j < n_bath_particles -1; ++j)
    {
        beta_a[j] = temp_a[j] - (g_dt * 0.25 * a_bath[j]) + (2 * g_dt * 0.25 * alpha_a[j]);
        beta_b[j] = temp_b[j] - (g_dt * 0.25 * b_bath[j]) + (2 * g_dt * 0.25 * alpha_b[j]);
    }

    // makes beta_a and beta_b the vector k_3
    f(beta_a,beta_b,current_time + (g_dt * 0.25));

    for (long long int j = 0; j < n_bath_particles -1; ++j)
    {
        a_bath[j] = temp_a[j] + (g_dt * 0.25 * ( (1.0/6.0)*a_bath[j] + (2.0/3.0)*alpha_a[j] + (1.0/6.0)*beta_a[j] ));
        b_bath[j] = temp_b[j] + (g_dt * 0.25 * ( (1.0/6.0)*b_bath[j] + (2.0/3.0)*alpha_b[j] + (1.0/6.0)*beta_b[j] ));
    }



    delete[] temp_a;
    delete[] temp_b;
    delete[] alpha_a;
    delete[] alpha_b;
    delete[] beta_a;
    delete[] beta_b;
    temp_a = nullptr;
    temp_b = nullptr;
    alpha_a = nullptr;
    alpha_b = nullptr;
    beta_a = nullptr;
    beta_b = nullptr;
} 

void update_bath_time_zero_quarter(double* a_0, double* b_0, double* a_bath, double* b_bath, double current_time)
{
    double* temp_a = new double[n_bath_particles -1];
    double* temp_b = new double[n_bath_particles -1];
    double* alpha_a = new double[n_bath_particles -1];
    double* alpha_b = new double[n_bath_particles -1];
    double* beta_a = new double[n_bath_particles -1];
    double* beta_b = new double[n_bath_particles -1];
    

    for (long long int j = 0; j < n_bath_particles -1; ++j)
    {
        temp_a[j] = a_bath[j];
        temp_b[j] = b_bath[j];
    }

    // Makes a_bath and b_bath the vector k_1
    f_time_zero(a_bath, b_bath, current_time);

    for (long long int j = 0; j < n_bath_particles -1; ++j)
    {
        alpha_a[j] = temp_a[j] + (0.5 * g_dt * 0.25 * a_bath[j]);
        alpha_b[j] = temp_b[j] + (0.5 * g_dt * 0.25 * b_bath[j]);
    }

    // makes alpha_a and alpha_b the vector k_2
    f(alpha_a, alpha_b, current_time + (0.5 * g_dt * 0.25));

    for (long long int j = 0; j < n_bath_particles -1; ++j)
    {
        beta_a[j] = temp_a[j] - (g_dt * 0.25 * a_bath[j]) + (2 * g_dt * 0.25 * alpha_a[j]);
        beta_b[j] = temp_b[j] - (g_dt * 0.25 * b_bath[j]) + (2 * g_dt * 0.25 * alpha_b[j]);
    }

    // makes beta_a and beta_b the vector k_3
    f(beta_a,beta_b,current_time + (g_dt * 0.25));

    for (long long int j = 0; j < n_bath_particles -1; ++j)
    {
        a_bath[j] = temp_a[j] + g_dt * 0.25 * ( (1.0/6.0)*a_bath[j] + (2.0/3.0)*alpha_a[j] + (1.0/6.0)*beta_a[j] );
        b_bath[j] = temp_b[j] + g_dt * 0.25 * ( (1.0/6.0)*b_bath[j] + (2.0/3.0)*alpha_b[j] + (1.0/6.0)*beta_b[j] );
    }



    delete[] temp_a;
    delete[] temp_b;
    delete[] alpha_a;
    delete[] alpha_b;
    delete[] beta_a;
    delete[] beta_b;
    temp_a = nullptr;
    temp_b = nullptr;
    alpha_a = nullptr;
    alpha_b = nullptr;
    beta_a = nullptr;
    beta_b = nullptr;
} 


// Compute the noise for either bath at a given time_step
double F(double* a_0, double* b_0, double* a_bath, double* b_bath, double* bath_init_pos, double* bath_init_vel, long long int current_fine_step)
{
    double output { 0 };
    #pragma omp parallel for reduction(+:output)
        for (long long int j = 0; j < n_bath_particles-1; ++j) {
            output += a_bath[j]*bath_init_vel[j+1] + b_bath[j]*bath_init_pos[j+1];
            //std::cout << output << std::endl;
        }
    output += a_0[current_fine_step]*bath_init_vel[0] + b_0[current_fine_step]*bath_init_pos[0];
    return output;
}



// Define a trapezoidal rule (no time scaling)
double trap(double *vel, double *theta, long long int current_step, double theta_norm)
{
    long long int max_steps { std::min(max_trap_steps, current_step - 1) };
    double output { 0 };
    //omp_set_num_threads(NUM_THREADS);
    #pragma omp parallel for reduction(+:output)
        for ( long long int j = 1; j < max_steps; ++j ) {
            //output += beta[current_step - j] * pos[j];
            output += theta[j] * vel[current_step - 1 - j];
        }
    output += ( 0.5 * theta[0] * vel[current_step - 1] ) + ( 0.5 * theta[max_steps] * vel[current_step - 1 - max_steps ] );
    /*output += - ( 3 * theta[max_steps-1] * vel[current_step - max_steps + 1] 
                - 4 * theta[max_steps-2] * vel[current_step - max_steps + 2]
                +     theta[max_steps-3] * vel[current_step - max_steps + 3]
                )/24.0
              - ( 3 * theta[0] * vel[current_step - 1]
                - 4 * theta[1] * vel[current_step - 2]
                +     theta[2] * vel[current_step - 3]
                )/24.0 ; */
    output = output/theta_norm;
    /*if (current_step % 1000 == 0)
    {
        std::cout << std::setprecision(std::numeric_limits<double>::max_digits10);
        std::cout << current_step << std::endl;
        std::cout << output << std::endl;
    }*/
    return output;
}

void get_interp_coef(double *vel, long long int current_step, double *coeff)
{
    // Just using 3 interpolation points for now
    coeff[0] = vel[current_step -1];
    coeff[1] = (vel[current_step -1] - vel[current_step -2])/g_dt;
    double temp_dd_coef {(vel[current_step -2] - vel[current_step -3])/g_dt};
    coeff[2] = (coeff[1] - temp_dd_coef)/(2*g_dt);
    
}

// Define a trapezoidal rule for the half step update of RK3
double trap_xtra(double *vel, double *theta,long long int current_step, double theta_norm, double* coeff, double xtra_step_size, double *shifted_theta)
{
    
    double output{0};
    long long int max_steps { std::min(max_trap_steps, current_step - 1) };
    if (xtra_step_size == 0.5)
    {
        #pragma omp parallel for reduction(+:output)
            for ( long long int j = 1; j < max_steps; ++j ) {
                output += vel[current_step - 1 - j] * shifted_theta[j];
                //std::cout << theta_f(g_dt*(j + xtra_step_size)) - shifted_theta[j] << std::endl;
        }
        // Above goes to max steps since it really inst the end of the trap rule. Below uses theta(0) = 1
        // to avoid the singularity in the function definition.
        // estimate vel(t + h/2) using quadratic interpolation
        double xtra_vel { coeff[0] + coeff[1]*g_dt *(xtra_step_size) + coeff[2]*(xtra_step_size)*(xtra_step_size + 1)*pow(g_dt,2) };
        output += ( (0.5 * xtra_vel) + (0.5 * shifted_theta[0] * vel[current_step - 1] ) ) * xtra_step_size;
        output += (0.5 * vel[current_step -1 - max_steps]*theta_f(g_dt*(max_steps + xtra_step_size))); 
        output += (0.5 * shifted_theta[0] * vel[current_step - 1] );

    }
    else
    {
        #pragma omp parallel for reduction(+:output)
            for ( long long int j = 1; j < max_steps; ++j ) {
                output += vel[current_step - 1 - j] * theta[j+1];
                //std::cout << theta_f(g_dt*(j + xtra_step_size)) - theta[j+1] << std::endl;
            }
        // Above goes to max steps since it really inst the end of the trap rule. Below uses theta(0) = 1
        // to avoid the singularity in the function definition.
        // estimate vel(t + h/2) using quadratic interpolation
        double xtra_vel { coeff[0] + coeff[1]*g_dt *(xtra_step_size) + coeff[2]*(xtra_step_size)*(xtra_step_size + 1)*pow(g_dt,2) };
        output += ( (0.5 * xtra_vel) + (0.5 * theta[1] * vel[current_step - 1] ) ) * xtra_step_size;
        output += (0.5 * vel[current_step -1 - max_steps]*theta_f(g_dt*(max_steps + xtra_step_size))); 
        output += (0.5 * theta[1] * vel[current_step - 1] );

        // endpoint corrections
        // FIGURE OUT ENDPOINT CORRECTIONS IN FUTURE VERSION FOR NOW I WILL TAKE O(dt^2) error.
    //}
    }
    //omp_set_num_threads(NUM_THREADS);
    
    
    return output/theta_norm;

}



// Derivative of the potential function
double fpot(double x)
{
    //return g_k_sqrd*x;
    return -12*pow(x,-13) + 6*pow(x,-7) + 4*pow(x - a0,3); // last term added for particle confinement
}


void f_system(double *pos, double *first_vel, double *last_vel, double* int_velocities, long long int current_step, double *theta, double left_noise, double right_noise, double theta_norm)
{

    double* temp_pos = new double[n_particles];
    double temp_first_vel {0};
    double temp_last_vel {0};
    double* temp_int_velocities = new double[n_int_particles];
    
    temp_pos[0] = first_vel[current_step-1];
    for (int j = 0; j < n_particles - 2; ++j)
    {
        temp_pos[j+1] = int_velocities[j];
    }
    temp_pos[n_particles - 1] = last_vel[current_step -1];


    temp_first_vel = fpot(pos[1] - pos[0] + a0)  + ( (g_k_sqrd * left_noise) / sqrt(theta_norm) ) -g_k_sqrd*trap(first_vel,theta, current_step, theta_norm)*g_dt ;
    for (int j = 1; j < n_int_particles+1; ++j)
    {
        temp_int_velocities[j-1] = fpot(pos[j+1] - pos[j] + a0) -fpot(pos[j] - pos[j -1] + a0);
    }

    temp_last_vel = -fpot(pos[n_particles -1] - pos[n_particles -2] + a0) + ( (g_k_sqrd * right_noise) / sqrt(theta_norm) ) -g_k_sqrd*trap(last_vel,theta, current_step, theta_norm)*g_dt ;

    for (int j = 0; j < n_particles; ++j)
    {
        pos[j] = temp_pos[j];
    }
    first_vel[current_step -1] = temp_first_vel;
    last_vel[current_step -1] = temp_last_vel;
    for (int j = 0; j < n_int_particles; ++j)
    {
        int_velocities[j] = temp_int_velocities[j];
    }

    delete[] temp_pos;
    temp_pos = nullptr;
    delete[] temp_int_velocities;
    temp_int_velocities = nullptr;

}

void f_system_finer_grids(double *pos, double *first_vel, double *last_vel, double* int_velocities, long long int current_step, double *theta, double left_noise, double right_noise, double theta_norm, double alter_step_size)
{
    //current_step +=1;
    double* temp_pos = new double[n_particles];
    double temp_first_vel {0};
    double temp_last_vel {0};
    double* temp_int_velocities = new double[n_int_particles];
    
    temp_pos[0] = first_vel[current_step-1];
    for (int j = 0; j < n_particles - 2; ++j)
    {
        temp_pos[j+1] = int_velocities[j];
    }
    temp_pos[n_particles - 1] = last_vel[current_step -1];

    temp_first_vel = fpot(pos[1] - pos[0] + a0) + ( (g_k_sqrd * left_noise) / sqrt(theta_norm) ) -g_k_sqrd*trap(first_vel,theta, current_step, theta_norm)*g_dt*alter_step_size ;
    for (int j = 1; j < n_int_particles+1; ++j)
    {
        temp_int_velocities[j-1] = fpot(pos[j+1] - pos[j] + a0) -fpot(pos[j] - pos[j -1] + a0);
    }
    temp_last_vel = -fpot(pos[n_particles -1] - pos[n_particles -2] + a0) + ( (g_k_sqrd * right_noise) / sqrt(theta_norm) ) -g_k_sqrd*trap(last_vel,theta, current_step, theta_norm)*g_dt*alter_step_size ;

    for (int j = 0; j < n_particles; ++j)
    {
        pos[j] = temp_pos[j];
    }
    first_vel[current_step -1] = temp_first_vel;
    last_vel[current_step -1] = temp_last_vel;
    for (int j = 0; j < n_int_particles; ++j)
    {
        int_velocities[j] = temp_int_velocities[j];
    }


    delete[] temp_pos;
    temp_pos = nullptr;
    delete[] temp_int_velocities;
    temp_int_velocities = nullptr;

}




void f_system_xtra(double *pos, double *first_vel, double *last_vel, double* int_velocities, long long int current_step, 
                    double *theta, double left_noise, double right_noise, double theta_norm, double xtra_step_size, 
                    double *extern_temp_first_vel, double *extern_temp_last_vel,
                    double hold_first_vel, double hold_last_vel, double *shifted_theta)
{
    double* temp_pos = new double[n_particles];
    double temp_first_vel {0};
    double temp_last_vel {0};
    double* temp_int_velocities = new double[n_int_particles];
    int n_interpolation_points {3};
    double* coeff = new double[n_interpolation_points];
    
    temp_pos[0] = first_vel[current_step-1];
    for (int j = 0; j < n_particles - 2; ++j)
    {
        temp_pos[j+1] = int_velocities[j];
    }
    temp_pos[n_particles - 1] = last_vel[current_step -1];


    double hold_extern_temp_first_vel {extern_temp_first_vel[current_step - 1]};
    extern_temp_first_vel[current_step - 1] = hold_first_vel;
    get_interp_coef( extern_temp_first_vel, current_step, coeff);
    temp_first_vel = fpot(pos[1] - pos[0] + a0) + ( (g_k_sqrd * left_noise) / sqrt(theta_norm) ) -g_k_sqrd*trap_xtra(extern_temp_first_vel,theta, current_step, theta_norm, coeff, xtra_step_size, shifted_theta)*g_dt ;
    extern_temp_first_vel[current_step -1] = hold_extern_temp_first_vel;
    
    for (int j = 1; j < n_int_particles+1; ++j)
    {
        temp_int_velocities[j-1] = fpot(pos[j+1] - pos[j] + a0) -fpot(pos[j] - pos[j -1] + a0);
    }

    double hold_extern_temp_last_vel {extern_temp_last_vel[current_step - 1]};
    extern_temp_last_vel[current_step - 1] = hold_last_vel;
    get_interp_coef( extern_temp_last_vel, current_step, coeff);
    temp_last_vel = -fpot(pos[n_particles -1] - pos[n_particles -2] + a0) + ( (g_k_sqrd * right_noise) / sqrt(theta_norm) ) - g_k_sqrd*trap_xtra(extern_temp_last_vel,theta, current_step, theta_norm, coeff, xtra_step_size, shifted_theta)*g_dt ;
    extern_temp_last_vel[current_step - 1] = hold_extern_temp_last_vel;


    for (int j = 0; j < n_particles; ++j)
    {
        pos[j] = temp_pos[j];
        //std::cout << pos[j] << ' ' <<  temp_pos[j] << std::endl;
    }
    
    first_vel[current_step -1] = temp_first_vel;
    last_vel[current_step -1] = temp_last_vel;
    for (int j = 0; j < n_int_particles; ++j)
    {
        int_velocities[j] = temp_int_velocities[j];
        //int_velocities[j] = 0;
    }

    delete[] temp_pos;
    temp_pos = nullptr;
    delete[] temp_int_velocities;
    temp_int_velocities = nullptr;
    delete[] coeff;    
    coeff = nullptr;
}

void RK3_system(double *pos, double *first_vel, double *last_vel, double* int_velocities, long long int current_step, 
double *theta, double &left_noise, double &right_noise, double theta_norm, double *a_0, double *b_0, double *a_bath,
double *b_bath, double *l_bath_init_pos, double *l_bath_init_vel, double *r_bath_init_pos,
double *r_bath_init_vel, double *shifted_theta)
{
    current_step +=1;
    
    double* temp_pos = new double[n_particles];
    double* temp_int_vel = new double[n_int_particles];
    double* alpha_pos = new double[n_particles];
    double* first_alpha_vel = new double[current_step];
    double* alpha_vel = new double[n_int_particles];
    double* last_alpha_vel = new double[current_step];
    double* beta_pos = new double[n_particles];
    double* first_beta_vel = new double[current_step];
    double* beta_vel = new double[n_int_particles];
    double* last_beta_vel = new double[current_step];
    double* coeff = new double[3];
    
    
    temp_pos[0] = pos[0];
    for (int j = 1; j < n_particles - 1; ++j)
    {
        temp_pos[j] = pos[j];
        temp_int_vel[j-1] = int_velocities[j-1];
    }
    temp_pos[n_particles -1] = pos[n_particles -1];


    double hold_first_vel {first_vel[current_step - 1]};
    double hold_last_vel {last_vel[current_step - 1]};

    

    // This assumes the noise is already at the correct time step
    f_system( pos, first_vel, last_vel, int_velocities, current_step, theta, left_noise, right_noise, theta_norm);
    
    alpha_pos[0] = temp_pos[0] + 0.5 * g_dt * pos[0];
    for (int j = 1; j < n_particles -1; ++j)
    {
        alpha_pos[j] = temp_pos[j] + 0.5 * g_dt * pos[j];
        alpha_vel[j-1] = temp_int_vel[j-1] + 0.5 * g_dt * int_velocities[j-1];
    }
    alpha_pos[n_particles -1] = temp_pos[n_particles -1] + 0.5 * g_dt * pos[n_particles - 1];

    /*
    for (long long int j = 0; j < current_step -1; ++j)
    {
        first_alpha_vel[j] = temp_first_vel[j] ;//+ 0.5 * g_dt * first_vel[j];
        last_alpha_vel[j] = temp_last_vel[j] ;//+ 0.5 * g_dt * last_vel[j];
    }
    */
    first_alpha_vel[current_step -1] = hold_first_vel + 0.5 * g_dt * first_vel[current_step -1];
    last_alpha_vel[current_step - 1] = hold_last_vel + 0.5 * g_dt * last_vel[current_step - 1];

    // Advance the noise (determines noise at half step))
    double current_time {g_dt * (current_step - 1)};
    update_bath(a_0, b_0, a_bath, b_bath, current_time); 
    left_noise  = F( a_0,  b_0, a_bath, b_bath, l_bath_init_pos, l_bath_init_vel, 2*current_step - 1);
    right_noise  = F( a_0,  b_0, a_bath, b_bath, r_bath_init_pos, r_bath_init_vel, 2*current_step - 1);
    // Perform the half update

    // The half is for the added half step of the trapezoidal rule
    double xtra_step_size {0.5};
    f_system_xtra( alpha_pos, first_alpha_vel, last_alpha_vel, alpha_vel, current_step, theta, left_noise, right_noise, theta_norm, xtra_step_size, first_vel, last_vel, hold_first_vel, hold_last_vel, shifted_theta);

    current_time += 0.5 * g_dt;

    beta_pos[0] = temp_pos[0] - ( g_dt * pos[0] ) + (2 * g_dt * alpha_pos[0] );
    for (int j = 1; j < n_particles -1; ++j)
    {
        beta_pos[j] = temp_pos[j] - ( g_dt * pos[j] ) + (2 * g_dt * alpha_pos[j] );
        beta_vel[j-1] = temp_int_vel[j-1] - ( g_dt * int_velocities[j-1]) + (2 * g_dt * alpha_vel[j-1]); 
    }
    beta_pos[n_particles -1] = temp_pos[n_particles -1] - ( g_dt * pos[n_particles -1] ) + (2 * g_dt * alpha_pos[n_particles -1] );

    /*
    for (long long int j = 0; j < current_step -1; ++j)
    {
        first_beta_vel[j] = temp_first_vel[j] ;//- g_dt * first_vel[j] + 2 * g_dt * first_alpha_vel[j];
        last_beta_vel[j] = temp_last_vel[j] ;//- g_dt * last_vel[j] + 2 * g_dt * last_alpha_vel[j];
    }
    */
    first_beta_vel[current_step -1] = hold_first_vel - g_dt * first_vel[current_step -1] + 2 * g_dt * first_alpha_vel[current_step -1];
    last_beta_vel[current_step -1] = hold_last_vel - g_dt * last_vel[current_step -1] + 2 * g_dt * last_alpha_vel[current_step -1];


    // Advance the noise (determines noise last stage of RK3)
    update_bath(a_0, b_0, a_bath, b_bath, current_time); 
    left_noise  = F( a_0,  b_0, a_bath, b_bath, l_bath_init_pos, l_bath_init_vel, 2*current_step);
    right_noise  = F( a_0,  b_0, a_bath, b_bath, r_bath_init_pos, r_bath_init_vel, 2*current_step);

    xtra_step_size = 1.0;

    f_system_xtra( beta_pos, first_beta_vel, last_beta_vel, beta_vel, current_step, theta, left_noise, right_noise, theta_norm, xtra_step_size, first_vel, last_vel, hold_first_vel, hold_last_vel, shifted_theta);


    pos[0] = temp_pos[0] + g_dt * ( (1.0/6.0) * pos[0] + (2.0/3.0) * alpha_pos[0] + (1.0/6.0) * beta_pos[0] );
    for (int j = 1; j < n_particles -1; ++j)
    {
        pos[j] = temp_pos[j] + g_dt * ( (1.0/6.0) * pos[j] + (2.0/3.0) * alpha_pos[j] + (1.0/6.0) * beta_pos[j] );
        int_velocities[j-1] = temp_int_vel[j-1] + g_dt * ( (1.0/6.0) * int_velocities[j-1] + (2.0/3.0) * alpha_vel[j-1] + (1.0/6.0) * beta_vel[j-1] );
    }
    pos[n_particles -1] = temp_pos[n_particles-1] + g_dt * ( (1.0/6.0) * pos[n_particles-1] + (2.0/3.0) * alpha_pos[n_particles-1] + (1.0/6.0) * beta_pos[n_particles-1] );

    

    // Possible overflow here?
    first_vel[current_step] = hold_first_vel + g_dt * ( (1.0/6.0) * first_vel[current_step -1] + (2.0/3.0) * first_alpha_vel[current_step -1] + (1.0/6.0) * first_beta_vel[current_step -1] );
    last_vel[current_step] = hold_last_vel + g_dt * ( (1.0/6.0) * last_vel[current_step -1] + (2.0/3.0) * last_alpha_vel[current_step -1] + (1.0/6.0) * last_beta_vel[current_step -1] );
    
    first_vel[current_step - 1]= hold_first_vel;
    last_vel[current_step - 1] = hold_last_vel;

    
    delete[] temp_pos;
    delete[] temp_int_vel;
    delete[] alpha_pos;
    delete[] first_alpha_vel;
    delete[] alpha_vel;
    delete[] last_alpha_vel;
    delete[] beta_pos;
    delete[] first_beta_vel;
    delete[] beta_vel;
    delete[] last_beta_vel;
    delete[] coeff;
    temp_pos = nullptr;
    temp_int_vel = nullptr;
    alpha_pos = nullptr;
    first_alpha_vel = nullptr;
    alpha_vel = nullptr;
    last_alpha_vel = nullptr;
    beta_pos = nullptr;
    first_beta_vel = nullptr;
    beta_vel = nullptr;
    last_beta_vel = nullptr;
    coeff = nullptr;
    
}

