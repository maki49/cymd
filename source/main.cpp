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
    //geo.write_coord("geo.out");

    //week 3,6: test search_adj
    /*
    clock_t start, end;
    start = clock();
    geo_init.search_adj(
        input.rcut_potential + input.extra_rcut_neighbor,
        input.max_neighbor);
    geo_init.search_adj_faster(
        input.rcut_potential + input.extra_rcut_neighbor,
        input.max_neighbor);
    end = clock();
    std::cout<<"time to search:  "<<double(end-start)/CLOCKS_PER_SEC<<"s"<<std::endl;  //unit: s
    geo_init.print_adj_list(12);
    */
    
    //week 7: iteration
    std::cout << input.nstep << std::endl;
    MD_step MS(input);
    MS.main_step(geo_init, lj_init);

    return 0;
}