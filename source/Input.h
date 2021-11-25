#ifndef _INPUT_
#define _INPUT_
#include <string>
class Input
{
public:
    Input();
    ~Input();
    void read_input_file(std::string& input_file);
    std::string geo_in_file;
    int natom;
    double rcut_potential;    //rcut for Leonard-Jones potential
    double extra_rcut_neighbor;   //extra distance for neighbor list 
    int max_neighbor;
    int whose_adj;  //which atom's adjacent atoms will be printed

    //parameters for L-J potential
    double epsilon;
    double sigma;

    double mass;
    int read_vel = 0;
    double init_temperature;  //(K)
    bool append = true;

    int nstep = 0;
    double dt = 0.01;//ps
    int steps_per_search = 5;   //howmany steps between each update of adjacent list
    int steps_per_print = 2;   // howmany steps between each print
    std::string ensemble = "NVE";
    int verlet_method = 1;
    
    //NVT paras
    double thermo_temperature;
    double nraise;

};
#endif
