#pragma once

class Noise
{
private:
    double* m_samples;
    double* m_s_vector;
    long long int m_samples_length;
    long long int m_s_vector_length;
    double m_temp;
    int m_t_f;
    double m_dt;

public:
    Noise( const double temp, const int t_f, const double dt);
    ~Noise();

    double& operator[] (const long long int current_step) const;
    long long int length() const;

};