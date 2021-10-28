#ifndef _LJ_POT_
#define _LJ_POT_
#include <vector>
#include <cmath>
#include "Geo.h"

class LJ_pot
{

public:
    double sigma;
    double epsilon;

    double Ep;
    std::vector<double> atom_pot;
    std::vector<vec3> atom_force;
    double ene_shift;   //u(rcut)
    double rcut;    //rcut_potential

    LJ_pot();
    LJ_pot(double sigma, double epsilon, double rcut);
    ~LJ_pot();

    void cal_EpF(Geo geo);
    void print_EF(int nat, int precision);
    vec3 f_ij(vec3 r12, double rcut);

private:
    double v_ij(vec3 r12, double rcut);

    double V_at(int ia, std::vector<int> adj_list_i, std::vector<vec3> adj_dis_list_i, double rcut);
    vec3 F_at(int ia, std::vector<int> adj_list_i, std::vector<vec3> adj_dis_list_i, double rcut);

};
#endif