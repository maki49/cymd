#ifndef _PS_
#define _PS_
#include "./Geo.h"
#include "./LJ_pot.h"
class Print_step
{
public:
    Print_step();
    Print_step(Geo& geo, LJ_pot& lj);
    ~Print_step();
    void print_info(int istep);
private:
    void print_position(int istep);
    void print_velocity(int istep);
    void print_force(int istep);
    void print_log(int istep);
    Geo geo_step;
    LJ_pot lj_step;
};
#endif