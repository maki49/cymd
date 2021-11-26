#ifndef _MD_
#define _MD_

#include "Geo.h"
#include "LJ_pot.h"
#include "Print_step.h"
#include "Input.h"
#include <random>

class MD_step
{
public:
    MD_step();
    MD_step(Input& input);
    ~MD_step();

    void main_step(Geo& geo_step, LJ_pot& lj_step);
    
private:
    double dt;
    double m;
    int nstep=0;
    int steps_per_print;
    int steps_per_search;
    double rcut_neighbor;
    //rcut_potential is fixed in lj_step;
    int max_neighbor;
    int verlet_method = 1;
    bool append = true;
    std::string ensemble;
    double thermo_temperature=0.0;
    double nraise;
    bool cal_msd=false;
    int msd_print_interval;

    //temp r, v, f (maybe change these to pointers!)
    //real memory
    std::vector<vec3> r1;
    std::vector<vec3> r2;
    std::vector<vec3> atom_msd;
    //pointers
    std::vector<vec3>* r_tmdt = &r1;
    std::vector<vec3>* r_tpdt = &r2;
    
    void allocate(int natom);
    void verlet_0(int istep, Geo& geo_step, LJ_pot& lj_step);
    void verlet_1(int istep, Geo& geo_step, LJ_pot& lj_step);
    void velocity_verlet_before(Geo& geo_step, LJ_pot& lj_step);
    void velocity_verlet_after(Geo& geo_step, LJ_pot& lj_step);
    void Anderson(Geo& geo_step, double sgm, 
        std::default_random_engine& generator);
    void print_msd(int istep, int natom, int precision);
};
#endif