#ifndef _MD_
#define _MD_

#include "Geo.h"
#include "LJ_pot.h"
#include "Print_step.h"
#include "Input.h"
class MD_step
{
public:
    MD_step();
    MD_step(Input& input);
    ~MD_step();
    
    void verlet(int istep, Geo& geo_step, LJ_pot& lj_step);
    void velocity_verlet(int istep, Geo& geo_step, LJ_pot& lj_step);
    void main_step(Geo& geo_step, LJ_pot& lj_step);
private:
    double dt;
    
    int nstep=0;
    int steps_per_print;
    int steps_per_search;
    double rcut_neighbor;
    //rcut_potential is fixed in lj_step;
    int max_neighbor;

    //temp r, v, f
    std::vector<vec3> r_tmdt;
    std::vector<vec3> r_tpdt;
    std::vector<vec3> r_tpdt_nobox; //used to calculate v

};
#endif