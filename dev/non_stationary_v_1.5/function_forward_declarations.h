#pragma once

// Forward declarations of the functions

double trap(double *vel, double *theta, long long int current_step, double theta_norm);
double trap_xtra_half(double *vel, double *theta,long long int current_step, double theta_norm);
double fpot(double x);
double beta_f(double t);
double a_0_f(double t);
double b_0_f(double t);
void update_bath(double* a_0, double* b_0, double* a_bath, double* b_bath, double current_time);
void update_bath_time_zero(double* a_0, double* b_0, double* a_bath, double* b_bath, double current_time);
double F(double* a_0, double* b_0, double* a_bath, double* b_bath, double* bath_init_pos, double* bath_init_vel, long long int current_fine_step);
double normalize_theta(double* theta, double scale_dt);
double theta_f(double t);


void RK3_system(double *pos, double *first_vel, double *last_vel, double* int_velocities, long long int current_step, 
    double *theta, double &left_noise, double &right_noise, double theta_norm, double *a_0, double *b_0, double *a_bath,
    double *b_bath, double *l_bath_init_pos, double *l_bath_init_vel, double *r_bath_init_pos,
    double *r_bath_init_vel, double *shifted_theta);

void f_system_finer_grids(double *pos, double *first_vel, double *last_vel, double* int_velocities, long long int current_step, double *theta,          double left_noise, double right_noise, double theta_norm, double alter_step_size);

void update_bath_time_zero_quarter(double* a_0, double* b_0, double* a_bath, double* b_bath, double current_time);

void update_bath_quarter(double* a_0, double* b_0, double* a_bath, double* b_bath, double current_time);