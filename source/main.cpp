#include <string>
#include <vector>
#include <fstream>
#include <iostream>
#include "./Geo.h"
#include "./Input.h"
#include "./LJ_pot.h"
#include "./Print_step.h"
#include "./MD_step.h"

int main(int argc, char* argv[])
{
    std::string input_file(argv[1]);

    Input input;
    input.read_input_file(input_file);
    
    Geo geo_init(input.mass, input.init_temperature);
    geo_init.read_geo(input.geo_in_file, input.read_vel);  //read coord and velocity
    
    LJ_pot lj_init(geo_init.natom,
        input.sigma, input.epsilon, input.rcut_potential);
    
    //week 7, 8: iteration
    MD_step MS(input);
    MS.main_step(geo_init, lj_init);

    return 0;
}