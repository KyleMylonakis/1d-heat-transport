#include <iostream>
#include "Weights_and_Nodes.h"
#include <algorithm>
#include <iomanip>

double fast_conv(Weights_and_Nodes &obj, double* vel, long long int current_step)
{

}


int main(int argc, char const *argv[])
{
   std::cout << std::setprecision(std::numeric_limits<double>::max_digits10);
    Weights_and_Nodes gamma_R {53};
    //gamma_R.print_nodes();
    //gamma_R.print_weights();
    //gamma_R.print_c();

    const long long int N {100*10000 + 1};
    const double dt {0.01};
    double* vel = new double[N];
    for( long long int jj = 0; jj < N; ++jj)
    {
        vel[jj] = sin(jj*dt);
        //vel[jj] = 1.0;
    }


    for( long long int jj = 0; jj < N-1; ++jj)
    {
        //std::cout << ones[jj] << std::endl;
        //std::cout << gamma_R.fast_conv() << std::endl;
        if( jj % 100 == 0)
        {
            std::cout << jj/100 << std::endl;
        }
        gamma_R.update_c(vel, jj + 1, dt);
        //std::cout << gamma_R.soe_approx(jj*dt) << std::endl;
    }
    std::cout << gamma_R.fast_conv() << std::endl;

    delete[] vel;
    vel = nullptr;
    return 0;
}
