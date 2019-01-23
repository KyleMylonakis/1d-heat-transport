#include "Noise.h"
#include "function_forward_declarations.h"
#include <iostream>
//#include <complex.h>
#include <fftw3.h>
#include <random>


Noise::Noise(const double temp, const int t_f, const double dt, Weights_and_Nodes &wghts_nds) : 
    m_temp {temp},
    m_t_f {t_f},
    m_dt {dt}
{
    m_samples_length = static_cast<int>(1/m_dt) * m_t_f + 1;
    //m_samples_length = 6;
    m_s_vector_length = 2*m_samples_length - 2;
    m_samples = new double[m_samples_length]();
    m_s_vector = new double[m_s_vector_length]();

    std::cout << "Initialzing s vector" << std::endl;
    m_s_vector[0] = m_temp; //Since theta(0) = 1
    for (long long int i = 1; i < m_samples_length; ++i)
    {
        m_s_vector[i] = m_temp*theta_f(i * m_dt, wghts_nds);
    }
    for (long long int i = 0; i < m_samples_length - 2; ++i)
    {
        m_s_vector[m_s_vector_length -1 - i] = m_temp * theta_f(m_dt * (i + 1), wghts_nds );
    }

    fftw_complex *in, *out;
    fftw_plan p;

    in = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * m_s_vector_length);
    out = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * m_s_vector_length);
    p = fftw_plan_dft_1d(m_s_vector_length, in, out, FFTW_BACKWARD, FFTW_ESTIMATE);

    for (long long int i = 0; i < m_s_vector_length; ++i)
    {
        in[i][0] = m_s_vector[i];
        in[i][1] = 0;
    }

    std::cout << "Performing the FFT" << std::endl;
    fftw_execute(p);

    std::cout << "Normalziing and sampling noise" << std::endl;
    for (long long int i = 0; i < m_s_vector_length; ++i)
    {
        if (out[i][0] < 0)
        {
            m_s_vector[i] = 0;
        }
        else
        {
            m_s_vector[i] = sqrt(out[i][0]/m_s_vector_length);
        }   
    }

    std::random_device rd;
    std::mt19937 mersenne(rd());
    //std::mt19937 mersenne(0); // fixed seed for testing purposes
    std::normal_distribution<double> distribution(0,1.0);


    /*
    std::cout << "Testing RNG" << std::endl;
    double* test_rng = new double[m_s_vector_length];
    double avg {0};
    double stdev {0};
    for (long long int i  = 0; i < m_s_vector_length; ++i)
    {
        test_rng[i] = distribution(mersenne);
        avg += test_rng[i];
        stdev += test_rng[i]*test_rng[i];
    }

    std::cout << "The empirical mean and variance are " << avg/m_s_vector_length << " and " << stdev/m_s_vector_length << std::endl;

    delete[] test_rng;
    test_rng = nullptr;
    */

    for (long long int i = 0; i < m_s_vector_length; ++i)
    {
        in[i][0] = distribution(mersenne)*m_s_vector[i];
        in[i][1] = distribution(mersenne)*m_s_vector[i];
    }

    std::cout << "Performing FFT" << std::endl;
    fftw_execute(p);


    for (long long int i = 0; i < m_samples_length; ++i)
    {
        m_samples[i] = out[i][0];
    }

    fftw_destroy_plan(p);
    fftw_free(in); 
    fftw_free(out);

    std::cout << "Noise object successfully constructed" << std::endl;
    /*
    for (long long int i = 0; i < m_s_vector_length; ++i)
    {
        std::cout << m_s_vector[i] << " " ;
    }
    */

}
    
Noise::~Noise()
{
    delete[] m_samples;
    delete[] m_s_vector;
    m_samples = nullptr;
    m_s_vector = nullptr;
    std::cout << "Noise Object Successfully Destructed" << std::endl;
}

double& Noise::operator[] (const long long int current_step) const
{
    return m_samples[current_step];
}

long long int Noise::length() const
{
    return m_samples_length;
}