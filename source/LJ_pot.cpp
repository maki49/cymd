#include <vector>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <assert.h>
#include "LJ_pot.h"

LJ_pot::LJ_pot() {}
LJ_pot::LJ_pot(int natom, double sigma, double epsilon, double rcut, vec3 R) :
    sigma(sigma), epsilon(epsilon), rcut(rcut), R(R)
{
    double sr6 = pow(this->sigma / this->rcut, 6);
    this->ene_shift = 2 * this->epsilon * sr6 * (sr6 - 1);
    this->atom_force.resize(natom);
    this->atom_pot.resize(natom);
}
LJ_pot::~LJ_pot(){}

double LJ_pot::v_ij(vec3 r12, double rcut)  //potential between 2 atoms
{
    if(r12.norm < rcut)
    {
        double sr6 = pow(this->sigma/r12.norm, 6);
        return 2*this->epsilon*sr6*(sr6-1)-this->ene_shift;
    }
    else
    {
        return 0;
    }
}


vec3 LJ_pot::f_ij(vec3 r12, double rcut)//force from atom j to atom i
{ 
    // 1=i, 2= j. Fij = -Fji
    if(r12.norm < rcut)
    {
        double sr6=pow(this->sigma/r12.norm, 6);
        return  r12*4*this->epsilon*6*sr6*
        (2*sr6-1)/r12.norm/r12.norm;   //r12 is (ri - rj), return Fi, so here is "+1"
    }
    else
    {
        return vec3(0, 0, 0);
    }
}

//calculate V at each atom's position 
double LJ_pot::V_at(int ia, std::vector<int> adj_list_i, std::vector<vec3>* atom_r, double rcut)
{
    double v = 0;
    for (int j =0;j<adj_list_i.size();++j)
    {
        v += this->v_ij(Geo::shortest(atom_r->at(ia)-atom_r->at(adj_list_i[j]), R), rcut);
    }
    return v;
}
//calculate F of each atom
vec3 LJ_pot::F_at(int ia, std::vector<int> adj_list_i, std::vector<vec3>* atom_r, double rcut)
{
    vec3 f(0,0,0);
    for (int j =0;j<adj_list_i.size();++j)
    {
        f += this->f_ij(Geo::shortest(atom_r->at(ia) - atom_r->at(adj_list_i[j]), R), rcut);
    }
    return f;
}
void LJ_pot::cal_EpF(Geo geo)
{
    vec3 sum_all_f(0, 0, 0);
    this->Ep=0;
    for (int ia = 0;ia<geo.natom;++ia)
    {
        this->atom_pot[ia]=V_at(ia, geo.adj_list[ia], geo.r_t, this->rcut);
        //this->Ep+=V_at(ia, geo.adj_list[ia], geo.adj_dis_list[ia], rcut);
        this->Ep+=this->atom_pot[ia];
        this->atom_force[ia]=F_at(ia,  geo.adj_list[ia], geo.r_t, this->rcut);
        sum_all_f += atom_force[ia];
    }
    assert(abs(sum_all_f.norm) < 1e-5);  //check force
    return;
}
void LJ_pot::print_EF(int nat, int precision)
{
    std::ofstream ofs;
    ofs.open("energy.txt");
    ofs <<"total energy (eV): "<<std::setprecision(precision)<<this->Ep<<std::endl;
    ofs<< "energy of each atom (eV): "<<std::endl;
    ofs<< std::setw(4)<<"num"<<std::setw(20)<<"energy(eV)"<<std::endl;
    for (int ia=0;ia<nat;++ia)
    {
        ofs<<std::setprecision(precision)<<std::setw(4)<<ia+1<<std::setw(20)<< this->atom_pot[ia]<<std::endl;
    }
    ofs.close();

    /*
    ofs.open("force.txt");
    ofs<< std::setw(4)<<"num"<<std::setw(20)<<"Fx(eV/Ang.)"<<std::setw(20)<<"Fy(eV/Ang.)"<<std::setw(20)<<"Fz(eV/Ang.)"<<std::endl;
    for (int ia=0;ia<nat;++ia)
    {
        ofs<<std::setprecision(precision)<<std::setw(4)<<ia+1<<std::setw(20)<< this->atom_force[ia].x
        <<std::setw(20)<<this->atom_force[ia].y<<std::setw(20)<<this->atom_force[ia].z<<std::endl;
    }
    ofs.close();
*/
}
