/*
* This file contains all the auxillary functions for the RK3 integrators
*/

//#include "constant_forward_declarations.h"

#include <iostream>
#include <cmath>
#include <algorithm>
#include <limits>
//#include <omp.h>
#include <iomanip>
#include "Noise.h"
#include "Weights_and_Nodes.h"

//#define NUM_THREADS 1

extern const double g_dt; // Time step
extern const double g_k_sqrd; // Coming from the second derivative of the Lennard Jones 
extern const double g_k;
extern const double a0; // Equilibrium spacing
extern const int n_particles;      // The number of particles in the resolved system
const int n_int_particles = n_particles - 2;
//extern const int n_bath_particles;
extern const long long int max_trap_steps;
extern const long long int total_steps;

extern const double r;
extern const double epsilon;
extern const double Tc;
extern const double T;


double cutoff_f(double t)
{
    return pow(r,(t/Tc)*(t/Tc));
}


// Define the function theta
double theta_f(double t, Weights_and_Nodes &wght_nds)
{
    //return cutoff_f(t) * (std::cyl_bessel_j(1,2*g_k*t))/(g_k*t);
    //return cutoff_f(t) * wght_nds.soe_approx(t);
    return wght_nds.soe_approx(t);
}

double normalize_theta(double *theta, double scale_dt)
{
    /*
    if (max_trap_steps >= total_steps )
    {
        return 1.0;
    }
    else
    {*/

        double output { 0 };
        long long int max_steps {std::min(max_trap_steps,total_steps)};
        for ( long long int j = 1; j < max_steps-1 ; ++j ) 
            {
                output += theta[j];
            }
        output += ( 0.5 * theta[0] ) + ( 0.5 * theta[max_steps - 1] );

        output += - ( 3 * theta[max_steps-1] //* vel[0] 
                - 4 * theta[max_steps-2] //* vel[1]
                +     theta[max_steps-3] //* vel[2]
                )/24.0
              - ( 3 * theta[0] //* vel[max_trap_steps - 1]
                - 4 * theta[1] //* vel[max_trap_steps - 2]
                +     theta[2] //* vel[max_trap_steps - 3]
                )/24.0 ; 

        output = scale_dt*g_dt*output;
        output = output*g_k;
        return output;
    //}
    
}

// Define a trapezoidal rule (no time scaling)
double trap(double *vel, double *theta, long long int current_step, double theta_norm, Weights_and_Nodes &wghts_nds)
{
    long long int max_steps {  current_step - 1 };
    double output { 0 };

    //std::cout << "The number of steps in the trap rule is " << max_steps << std::endl;

    /*
    // For debugging
    for (long long int i = 0; i < max_steps + 1; ++i)
    {
     vel[i] =pow(i*g_dt,2);
    }
    */

    /*
    if (current_step < 4)
    {
        #pragma omp parallel for reduction(+:output)
        for ( long long int j = 1; j < max_steps; ++j ) {
            //output += beta[current_step - j] * pos[j];
            output += theta[j] * vel[current_step - 1 - j];
        }
        output += ( 0.5 * theta[0] * vel[current_step - 1] ) + ( 0.5 * theta[max_steps] * vel[current_step - 1 - max_steps ] );
        return output/theta_norm;
    }
    //omp_set_num_threads(NUM_THREADS);
    /*
    #pragma omp parallel for reduction(+:output)
        for ( long long int j = 1; j < max_steps; ++j ) {
            //output += beta[current_step - j] * pos[j];
            output += theta[j] * vel[current_step - 1 - j];
        }
    output += ( 0.5 * theta[0] * vel[current_step - 1] ) + ( 0.5 * theta[max_steps] * vel[current_step - 1 - max_steps ] );
    */
    output = wghts_nds.fast_conv();

    /*
    output += - ( 3 * theta[max_steps] * vel[current_step - max_steps - 1] 
                - 4 * theta[max_steps - 1] * vel[current_step - max_steps]
                +     theta[max_steps - 2] * vel[current_step - max_steps + 1]
                )/24.0
              - ( 3 * theta[0] * vel[current_step - 1]
                - 4 * theta[1] * vel[current_step - 2]
                +     theta[2] * vel[current_step - 3]
                )/24.0 ;  
                */
    output = output/theta_norm;
    /*
    if ((current_step-1) % static_cast<int>(1/g_dt) == 0)
    {
        std::cout << std::setprecision(std::numeric_limits<double>::max_digits10);
        std::cout << "The endpoint corrected trapezloidal rule has value " << output*g_dt <<
                        " at time " << (current_step-1)*g_dt << std::endl;
    }
    */
    return output;
}

void get_interp_coef(double *vel, long long int current_step, double *coeff)
{
    // For debugging
    /*
    long long int max_steps { std::min(max_trap_steps, current_step - 1) };
    for (long long int i = 0; i < max_steps + 1; ++i)
    {
     vel[i] =pow(i*g_dt,2);
    }
    */


    // Just using 3 interpolation points for now
    coeff[0] = vel[current_step -1];
    coeff[1] = (vel[current_step -1] - vel[current_step -2])/g_dt;
    double temp_dd_coef {(vel[current_step -2] - vel[current_step -3])/g_dt};
    coeff[2] = (coeff[1] - temp_dd_coef)/(2*g_dt);
    
}

// Define a trapezoidal rule for the half step update of RK3
double trap_xtra(double *vel, double *theta,long long int current_step, double theta_norm, double* coeff, double xtra_step_size, double *shifted_theta,
                    Weights_and_Nodes &wghts_nds)
{
    
    //double output{0};
    long long int max_steps { current_step - 1 };
    double xtra_vel { coeff[0] + coeff[1]*g_dt *(xtra_step_size) + coeff[2]*(xtra_step_size)*(xtra_step_size + 1)*pow(g_dt,2) };
    double output {wghts_nds.twisted_conv(xtra_step_size * g_dt)};
    output += 0.5* (xtra_step_size*g_dt) * (wghts_nds.soe_approx(0)*xtra_vel + wghts_nds.soe_approx(xtra_step_size*g_dt)*vel[current_step - 1] );
    return output/theta_norm;

    // For debugging
    //max_steps = (((3*3*3*11*11*11)/3) + 1)/2;
    //current_step = (((3*3*3*11*11*11)/3) + 1)/2;
    //std::cout << "The current time is " << (current_step -1)*g_dt << std::endl;
    /*
    for (long long int i = 0; i < max_steps + 1; ++i)
    {
     vel[i] =pow(i*g_dt,2);
    }

   /*
    if (xtra_step_size == 0.5)
    {
        double xtra_vel { coeff[0] + coeff[1]*g_dt *(xtra_step_size) + coeff[2]*(xtra_step_size)*(xtra_step_size + 1)*pow(g_dt,2) };
        double output {wghts_nds.soe_approx(0.5 * dt)* wghts_nds.fast_conv()};
        output += 0.5* (0.5*dt) * (wghts_nds.soe_approx(0)*xtra_vel + wghts_nds.soe_approx(0.5*dt)*vel[current_step - 1] );
        return output/trunc_theta_norm;
        /*
        #pragma omp parallel for reduction(+:output)
            for ( long long int j = 1; j < max_steps; ++j ) {
                output += vel[current_step - 1 - j] * shifted_theta[j];
                //std::cout << theta_f(g_dt*(j + xtra_step_size)) - shifted_theta[j] << std::endl;
        }
        // Above goes to max steps since it really inst the end of the trap rule. Below uses theta(0) = 1
        // to avoid the singularity in the function definition.
        // estimate vel(t + h/2) using quadratic interpolation
        double xtra_vel { coeff[0] + coeff[1]*g_dt *(xtra_step_size) + coeff[2]*(xtra_step_size)*(xtra_step_size + 1)*pow(g_dt,2) };
        double xtra_vel_minus { coeff[0] + coeff[1]*g_dt *(-xtra_step_size) + coeff[2]*(-xtra_step_size)*(-xtra_step_size + 1)*pow(g_dt,2) };

        output += ( (0.5 * xtra_vel) + (0.5 * shifted_theta[0] * vel[current_step - 1] ) ) * xtra_step_size;
        
        output += ( - ( -0.5 * theta[0]* xtra_vel + 0.5 * theta[1]*xtra_vel_minus )/12.0
              - ( 3 * xtra_vel*theta[0]
                - 4 * vel[current_step -1]*shifted_theta[0]
                +     xtra_vel_minus*theta[1]
                )/24.0 ) * xtra_step_size;


        output += (0.5 * vel[current_step -1 - max_steps]*theta_f(g_dt*(max_steps + xtra_step_size))); 
        output += (0.5 * shifted_theta[0] * vel[current_step - 1] );
        
        output += - ( 3 * vel[current_step -1 - max_steps]*theta_f(g_dt*(max_steps + xtra_step_size)) 
                - 4 * vel[current_step -1 - max_steps + 1]*theta_f(g_dt*(max_steps + xtra_step_size - 1))
                +     vel[current_step -1 - max_steps + 2]*theta_f(g_dt*(max_steps + xtra_step_size - 2))
                )/24.0
              - ( 3 * vel[current_step -1]*shifted_theta[0]
                - 4 * vel[current_step -2]*shifted_theta[1]
                +     vel[current_step -3]*shifted_theta[2]
                )/24.0 ;  
        /*
        if ( std::abs( (current_step-1 + xtra_step_size)*g_dt - 1.00 ) < pow(10,-12))
        {
            std::cout << std::setprecision(std::numeric_limits<double>::max_digits10);
            std::cout << "The endpoint corrected trapezloidal rule has value " << output*g_dt <<
                            " at time " << (current_step-1 + xtra_step_size)*g_dt << std::endl;
        }
        */
       
        
/*
    }
    else
    {
        double xtra_vel { coeff[0] + coeff[1]*g_dt *(xtra_step_size) + coeff[2]*(xtra_step_size)*(xtra_step_size + 1)*pow(g_dt,2) };
        double output {wghts_nds.soe_approx(xtra_step_size * dt)* wghts_nds.fast_conv()};
        output += 0.5* (xtra_step_size*dt) * (wghts_nds.soe_approx(0)*xtra_vel + wghts_nds.soe_approx(xtra_step_size*dt)*vel[current_step - 1] );
        return output/trunc_theta_norm;
        /*
        #pragma omp parallel for reduction(+:output)
            for ( long long int j = 1; j < max_steps; ++j ) {
                output += vel[current_step - j] * theta[j];
                //std::cout << theta_f(g_dt*(j + xtra_step_size)) - theta[j+1] << std::endl;
            }
        // Above goes to max steps since it really inst the end of the trap rule. Below uses theta(0) = 1
        // to avoid the singularity in the function definition.
        // estimate vel(t + h/2) using quadratic interpolation
        double xtra_vel { coeff[0] + coeff[1]*g_dt *(xtra_step_size) + coeff[2]*(xtra_step_size)*(xtra_step_size + 1)*pow(g_dt,2) };
        output += ( (0.5 * xtra_vel) );//+ (0.5 * theta[1] * vel[current_step - 1] ) ) * xtra_step_size;
        output += (0.5 * vel[current_step -1 - max_steps]*theta_f(g_dt*(max_steps + xtra_step_size))); 
        //output += (0.5 * theta[1] * vel[current_step - 1] );
        
        
        output += - ( 3 * vel[current_step -1 - max_steps]*theta_f(g_dt*(max_steps + xtra_step_size)) 
                - 4 * vel[current_step -1 - max_steps + 1]*theta_f(g_dt*(max_steps + xtra_step_size - 1))
                +     vel[current_step -1 - max_steps + 2]*theta_f(g_dt*(max_steps + xtra_step_size - 2))
                )/24.0
              - ( 3 * xtra_vel*theta[0]
                - 4 * vel[current_step -1]*theta[1]
                +     vel[current_step -2]*theta[2]
                )/24.0 ; 

        /*
        if ( std::abs( (current_step-1 + xtra_step_size)*g_dt - 1.00 ) < pow(10,-8))
        {
            std::cout << std::setprecision(std::numeric_limits<double>::max_digits10);
            std::cout << "The endpoint corrected trapezloidal rule has value " << output*g_dt <<
                            " at time " << (current_step-1 + xtra_step_size)*g_dt << std::endl;
        }
        */
    //}
    
    //}
    //omp_set_num_threads(NUM_THREADS);
    
    
    //return output/theta_norm;


}



// Derivative of the potential function
double fpot(double x)
{
    //return g_k_sqrd*x;
    return -12*pow(x,-13) + 6*pow(x,-7) + 4*pow(x - a0,3); // last term added for particle confinement
}


void f_system(double *pos, double *first_vel, double *last_vel, double* int_velocities, 
    long long int current_step, double *theta, double left_noise, double right_noise, double theta_norm,
    Weights_and_Nodes &left_wghts_nds, Weights_and_Nodes &right_wghts_nds)
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


    temp_first_vel = fpot(pos[1] - pos[0] + a0)  + ( (g_k_sqrd * left_noise) ) -g_k_sqrd*trap(first_vel,theta, current_step, theta_norm, left_wghts_nds) ;
    for (int j = 1; j < n_int_particles+1; ++j)
    {
        temp_int_velocities[j-1] = fpot(pos[j+1] - pos[j] + a0) -fpot(pos[j] - pos[j -1] + a0);
    }

    temp_last_vel = -fpot(pos[n_particles -1] - pos[n_particles -2] + a0) + ( (g_k_sqrd * right_noise) ) -g_k_sqrd*trap(last_vel,theta, current_step, theta_norm, right_wghts_nds) ;

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


double old_trap(double *vel, double *theta, long long int current_step, double theta_norm)
{
    long long int max_steps { std::min(max_trap_steps, current_step - 1) };
    double output { 0 };

    //std::cout << "The number of steps in the trap rule is " << max_steps << std::endl;

    /*
    // For debugging
    for (long long int i = 0; i < max_steps + 1; ++i)
    {
     vel[i] =pow(i*g_dt,2);
    }
    */

    if (current_step < 4)
    {
        #pragma omp parallel for reduction(+:output)
        for ( long long int j = 1; j < max_steps; ++j ) {
            //output += beta[current_step - j] * pos[j];
            output += theta[j] * vel[current_step - 1 - j];
        }
        output += ( 0.5 * theta[0] * vel[current_step - 1] ) + ( 0.5 * theta[max_steps] * vel[current_step - 1 - max_steps ] );
        return output/theta_norm;
    }
    //omp_set_num_threads(NUM_THREADS);
    #pragma omp parallel for reduction(+:output)
        for ( long long int j = 1; j < max_steps; ++j ) {
            //output += beta[current_step - j] * pos[j];
            output += theta[j] * vel[current_step - 1 - j];
        }
    output += ( 0.5 * theta[0] * vel[current_step - 1] ) + ( 0.5 * theta[max_steps] * vel[current_step - 1 - max_steps ] );
    
    output += - ( 3 * theta[max_steps] * vel[current_step - max_steps - 1] 
                - 4 * theta[max_steps - 1] * vel[current_step - max_steps]
                +     theta[max_steps - 2] * vel[current_step - max_steps + 1]
                )/24.0
              - ( 3 * theta[0] * vel[current_step - 1]
                - 4 * theta[1] * vel[current_step - 2]
                +     theta[2] * vel[current_step - 3]
                )/24.0 ;  
    output = output/theta_norm;
    /*
    if ((current_step-1) % static_cast<int>(1/g_dt) == 0)
    {
        std::cout << std::setprecision(std::numeric_limits<double>::max_digits10);
        std::cout << "The endpoint corrected trapezloidal rule has value " << output*g_dt <<
                        " at time " << (current_step-1)*g_dt << std::endl;
    }
    */
    return output;
}




void f_system_finer_grids(double *pos, double *first_vel, double *last_vel, double* int_velocities, long long int current_step, double *theta, double* left_noise, double* right_noise, double theta_norm, double alter_step_size)
{
    //current_step +=1;
    double* temp_pos = new double[n_particles];
    double temp_first_vel {0};
    double temp_last_vel {0};
    double* temp_int_velocities = new double[n_int_particles];

    /*
    temp_pos[0] = pos[0];
    pos[0] = -first_vel[current_step-1];
    first_vel[current_step - 1] = temp_pos[0];
    */
    
    
    temp_pos[0] = first_vel[current_step-1];
    for (int j = 0; j < n_particles - 2; ++j)
    {
        temp_pos[j+1] = int_velocities[j];
    }
    temp_pos[n_particles - 1] = last_vel[current_step -1];

    temp_first_vel = fpot(pos[1] - pos[0] + a0) + ( (g_k_sqrd * left_noise[current_step - 1]) ) -g_k_sqrd*old_trap(first_vel,theta, current_step, theta_norm)*g_dt*alter_step_size ;
    for (int j = 1; j < n_int_particles+1; ++j)
    {
        temp_int_velocities[j-1] = fpot(pos[j+1] - pos[j] + a0) -fpot(pos[j] - pos[j -1] + a0);
    }
    temp_last_vel = -fpot(pos[n_particles -1] - pos[n_particles -2] + a0) + ( (g_k_sqrd * right_noise[current_step - 1]) / sqrt(theta_norm) ) -g_k_sqrd*old_trap(last_vel,theta, current_step, theta_norm)*g_dt*alter_step_size ;

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
                    double hold_first_vel, double hold_last_vel, double *shifted_theta,
                    Weights_and_Nodes &left_wghts_nds, Weights_and_Nodes &right_wghts_nds)
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
    temp_first_vel = fpot(pos[1] - pos[0] + a0) + ( (g_k_sqrd * left_noise) ) -g_k_sqrd*trap_xtra(extern_temp_first_vel,theta, current_step, theta_norm, coeff, xtra_step_size, shifted_theta, left_wghts_nds) ;
    extern_temp_first_vel[current_step -1] = hold_extern_temp_first_vel;
    
    for (int j = 1; j < n_int_particles+1; ++j)
    {
        temp_int_velocities[j-1] = fpot(pos[j+1] - pos[j] + a0) -fpot(pos[j] - pos[j -1] + a0);
    }

    double hold_extern_temp_last_vel {extern_temp_last_vel[current_step - 1]};
    extern_temp_last_vel[current_step - 1] = hold_last_vel;
    get_interp_coef( extern_temp_last_vel, current_step, coeff);
    temp_last_vel = -fpot(pos[n_particles -1] - pos[n_particles -2] + a0) + ( (g_k_sqrd * right_noise)  ) - g_k_sqrd*trap_xtra(extern_temp_last_vel,theta, current_step, theta_norm, coeff, xtra_step_size, shifted_theta, right_wghts_nds) ;
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
    double *theta, double* left_noise, double* right_noise, double theta_norm, double *shifted_theta,
    Weights_and_Nodes &left_wghts_nds, Weights_and_Nodes &right_wghts_nds)
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
    
    left_wghts_nds.update_c(first_vel, current_step-1, g_dt);
    right_wghts_nds.update_c(last_vel, current_step-1, g_dt);
    
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
    f_system( pos, first_vel, last_vel, int_velocities, current_step, theta, left_noise[2*(current_step - 1)], right_noise[2*(current_step - 1)], theta_norm, left_wghts_nds, right_wghts_nds);
    
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
   
    // Perform the half update

    // The half is for the added half step of the trapezoidal rule
    double xtra_step_size {0.5};
    f_system_xtra( alpha_pos, first_alpha_vel, last_alpha_vel, alpha_vel, current_step, theta, left_noise[2*(current_step - 1) + 1], right_noise[2*(current_step - 1) + 1], 
                    theta_norm, xtra_step_size, first_vel, last_vel, hold_first_vel, hold_last_vel, shifted_theta, left_wghts_nds, right_wghts_nds);

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
    
    xtra_step_size = 1.0;

    f_system_xtra( beta_pos, first_beta_vel, last_beta_vel, beta_vel, current_step, theta, left_noise[2*(current_step)], right_noise[2*(current_step)], theta_norm, xtra_step_size, first_vel, last_vel, hold_first_vel, hold_last_vel, shifted_theta, left_wghts_nds, right_wghts_nds);


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

