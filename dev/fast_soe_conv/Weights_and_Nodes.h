#include <complex>
#include <fstream>


class Weights_and_Nodes
{
private:
    std::complex<double>* m_nodes;
    std::complex<double>* m_weights;
    int m_M;

    std::complex<double>* m_c;

public:

    Weights_and_Nodes(int M);
    ~Weights_and_Nodes();

    void print_nodes();
    void print_weights();
    void print_c();

    void update_c(double* vel, long long int current_step, double dt);
    double fast_conv();
    double twisted_conv(double dt);
    double soe_approx(double t);

};