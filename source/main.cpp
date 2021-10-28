#include <string>
#include <vector>
#include <fstream>
#include <iostream>
#include "./Geo.h"
#include "./Input.h"
#include "./LJ_pot.h"
#include "./Print_step.h"

int main(int argc, char* argv[])
{
    std::string input_file(argv[1]);
    //week 1, 5

    Input input;
    input.read_input_file(input_file);
    Geo geo(input.mass, input.init_temperature);
    
    geo.read_geo(input.geo_in_file, input.read_vel);  //read coord and velocity
    //std::cout << "read ok" << std::endl;
    //geo.write_coord("geo.out");
    clock_t start, end;
    //week 3
    start = clock();
    geo.search_adj(
        input.rcut_potential + input.extra_rcut_neighbor,
        input.max_neighbor);
    /*
    geo.search_adj_faster(
        input.rcut_potential + input.extra_rcut_neighbor,
        input.max_neighbor);
*/
    end = clock();
    std::cout<<"time to search:  "<<double(end-start)/CLOCKS_PER_SEC<<"s"<<std::endl;  //unit: s
    geo.print_adj_list(12);
/*
    //=======step 0===============
    //week 4
    LJ_pot lj(input.sigma, input.epsilon, input.rcut_potential);
    lj.cal_EpF(geo);
    //lj.print_EF(geo.natom, geo.precision);
    
    //week 5
    Print_step ps(geo, lj);
    ps.print_info(0, input.append);   //E, f, x, v of istep=0
    //=======\step 0===============
*/
    return 0;
}