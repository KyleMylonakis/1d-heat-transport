#include "Weights_and_Nodes.h"
#include <iostream>
#include <cmath>

Weights_and_Nodes::Weights_and_Nodes(const int M) : m_M {M}
{
    m_nodes = new std::complex<double>[m_M];
    m_weights = new std::complex<double>[m_M];
    m_c = new std::complex<double>[m_M]{};

    std::ifstream gamma_real_f;
    std::ifstream gamma_imag_f;
    std::ifstream R_real_f;
    std::ifstream R_imag_f;
    gamma_real_f.open("gamma_real.bin", std::ios::binary | std::ios::ate);
    gamma_imag_f.open("gamma_imag.bin", std::ios::binary | std::ios::ate);
    R_real_f.open("R_real.bin", std::ios::binary | std::ios::ate);
    R_imag_f.open("R_imag.bin", std::ios::binary | std::ios::ate);
    
    if (gamma_real_f.is_open() &&
        gamma_imag_f.is_open() &&
        R_real_f.is_open() &&
        R_imag_f.is_open())
        {
            std::streampos size_g_r { gamma_real_f.tellg() };
            std::streampos size_g_i { gamma_imag_f.tellg() };
            std::streampos size_r_r { R_real_f.tellg() };
            std::streampos size_r_i { R_imag_f.tellg() };

            char* memblock_g_r = new char[size_g_r];
            char* memblock_g_i = new char[size_g_i];
            char* memblock_r_r = new char[size_r_r];
            char* memblock_r_i = new char[size_r_i];

            gamma_real_f.seekg(0,std::ios::beg);
            gamma_imag_f.seekg(0,std::ios::beg);
            R_real_f.seekg(0,std::ios::beg);
            R_imag_f.seekg(0,std::ios::beg);

            gamma_real_f.read(memblock_g_r, size_g_r);
            gamma_imag_f.read(memblock_g_i, size_g_i);
            R_real_f.read(memblock_r_r, size_r_r);
            R_imag_f.read(memblock_r_i, size_r_i);

            gamma_real_f.close();
            gamma_imag_f.close();
            R_real_f.close();
            R_imag_f.close();

            // Reinterpret as doubles
            double* gamma_real_values = (double*)memblock_g_r;
            double* gamma_imag_values = (double*)memblock_g_i;
            double* R_real_values = (double*)memblock_r_r;
            double* R_imag_values = (double*)memblock_r_i; 
            
            
            
            for(int ii=0; ii<M; ++ii)
            {
            m_nodes[ii] = std::complex<double>( gamma_real_values[ii], gamma_imag_values[ii]);
            m_weights[ii] = { R_real_values[ii], R_imag_values[ii] };
            }

            std::cout << "Loaded weights and nodes of the sum of exponentials approximation." << std::endl;
        }
        else
        {
            std::cerr << "Failed to load weights and nodes" << std::endl;
        }
}

Weights_and_Nodes::~Weights_and_Nodes()
{
    delete[] m_weights;
    m_weights = nullptr;

    delete[] m_nodes;
    m_nodes = nullptr;

    std::cout << "Weights and Nodes Deallocated" << std::endl;
}

void Weights_and_Nodes::print_nodes()
{
    for( int jj = 0; jj < m_M; ++jj)
    {
        std::cout << "Node[" << jj + 1 << "]: " << m_nodes[jj] << std::endl;
    }
}

void Weights_and_Nodes::print_weights()
{
    for( int jj = 0; jj < m_M; ++jj)
    {
        std::cout << "weight[" << jj + 1 << "]: " << m_weights[jj] << std::endl;
    }
}

void Weights_and_Nodes::print_c()
{
    for( int jj = 0; jj < m_M; ++jj)
    {
        std::cout << "c[" << jj << "]: " << m_c[jj] << std::endl;
    }
}


void Weights_and_Nodes::update_c(double* vel, long long int current_step, double dt)
{
    if( current_step != 0)
    {
        std::complex<double> next_step {0};
        for( int jj = 0; jj < m_M; ++jj)
        {
            next_step = 0.5*dt*( vel[current_step] 
            + vel[current_step - 1]*std::exp( m_nodes[jj]*dt ) );
            m_c[jj] =  std::exp(m_nodes[jj]*dt)*m_c[jj] + next_step;
        }
    }
}




double Weights_and_Nodes::fast_conv()
{
    std::complex<double> out {0,0};
    for( int jj = 0; jj < m_M; ++jj)
    {
        out += m_weights[jj]*m_c[jj];
    }
    return out.real();
}

double Weights_and_Nodes::twisted_conv(double dt)
{
    std::complex<double> out {0,0};
    for (int jj =0; jj < m_M; ++jj)
    {
        out +=m_weights[jj]*std::exp(dt * m_nodes[jj])*m_c[jj];
    }
    return out.real();
}


double Weights_and_Nodes::soe_approx(double t)
{
    std::complex<double> out {0};
    for( int jj = 0; jj < m_M; ++jj)
    {
        out += m_weights[jj]*std::exp(m_nodes[jj]*t);
    }
    return out.real();
}