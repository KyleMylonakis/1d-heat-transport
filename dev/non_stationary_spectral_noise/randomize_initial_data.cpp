#include <random>
#include <iostream>
#include <ctime>

extern const double left_temp;
extern const double right_temp;
extern const double avg_init_temp;
extern const double a0;
extern const double max_rand_radius { std::min(avg_init_temp, a0/4.0)};

// Initialize Normal and Uniform RNG
std::random_device rd;
std::mt19937 generator(rd());
double mean { 0 };
double variance { avg_init_temp };
std::normal_distribution<double> norm_dist(mean, sqrt(variance)) ;
std::uniform_real_distribution<double> unif_dist(-max_rand_radius, max_rand_radius); 


void initialize_velocities(double* int_velocities, double* first_vel, double* last_vel, double num_int_particles)
{
    // Initializes velocities as i.i.d gaussians
    for (int ii = 0; ii < num_int_particles; ++ii)
    {
        int_velocities[ii] = norm_dist(generator);
    }
    first_vel[0] = norm_dist(generator);
    last_vel[0] = norm_dist(generator);
}

void initialize_positions(double* positions, int n_particles)
{
    // Initializes positions from unfirom distribution
    for (int ii = 0; ii < n_particles; ++ii)
    {
        positions[ii] = unif_dist(generator);
    }
}

void compute_mean(double* data, int num_samp)
{
    double mean {0};
    for (int ii = 0; ii < num_samp; ++ii)
    {
        mean += data[ii];
    }
    mean /= num_samp;
    std::cout << "The sample mean is " << mean << std::endl;
}

void compute_variance(double* data, int num_samp)
{
    double mean {0};
    for (int ii = 0; ii < num_samp; ++ii)
    {
        mean += data[ii]*data[ii];
    }
    mean /= num_samp;
    std::cout << "The sample variance is " << mean << std::endl;
}