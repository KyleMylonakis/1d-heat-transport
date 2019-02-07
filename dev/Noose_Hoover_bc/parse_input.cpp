/*
    * Function library responsible for converting 
    * the command line arguments to usable variables
*/

#include <iostream>
#include <string>
#include <sstream>
#include <cstdlib>

template <typename T>

void parse_inputs(char* argument, T& save_input)
{
    std::stringstream input_buffer { argument };
    if (!(input_buffer >> save_input))
    {
        std::cerr << "Error: Invalid input" << std::endl;
        exit(1);
    }
    std::cout << "Command line argument set as " << save_input << std::endl;
}

/*
void parse_inputs(char* argument, double& save_input)
{
    std::stringstream input_buffer { argument };
    if (!(input_buffer >> save_input))
    {
        std::cerr << "Error: Invalid input" << std::endl;
        exit(1);
    }
    std::cout << "Command line argument set as " << save_input << std::endl;
}

void parse_inputs(char* argument, int& save_input)
{
    std::stringstream input_buffer { argument };
    if (!(input_buffer >> save_input))
    {
        std::cerr << "Error: Invalid input" << std::endl;
        exit(1);
    }
    std::cout << "Command line argument set as " << save_input << std::endl;
}
*/
