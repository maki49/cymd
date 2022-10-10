#ifndef _GEO_
#define _GEO_
#include <fstream>
#include "vec3.h"
class Geo
{
public:
    const double kb = 1.38064852;    //*1e-23, bolzmann constant
    const double NA = 6.02214086;      //*1e23
    const double Ap2ms = 100;  //Ang/ps -> m/s
    const double e = 1.602176634;   //*1e-19 C

    Geo();
    Geo(double m, double T0);
    ~Geo();
    
    void read_geo(std::string& geo_in_file, int read_vel, std::string v0_type);
    //void write_coord(std::string& geo_out_file);
    //search adjacent atoms for adj_list, NOT including itself
    void search_adj(double rcut, int max_neighbor);
    //only for rcut < L/2, and L is equal in 3 dimensions
    void search_adj_faster(double rcut, int max_neighbor);
    void cal_EkT(double sum_v2);
    void print_adj_list(int ia);

    //called before calculating Ep and Force, 
    //update  according to the existing adj_list
    void update_dis_list(void);
    static vec3 res_in_box(vec3 r, vec3 R);
    static vec3 shortest(vec3 r, vec3 R);

    int natom = 0;
    vec3 lattice_vec[3];
    vec3 R;
    int precision = 0;  //decimal places
    std::vector<std::string> atom_types;
    std::vector<vec3> atom_coords;
    std::vector<vec3>* r_t = &atom_coords; //pointer at `atom_coords` at time t
    std::vector<vec3> atom_v;   //Ang./ps
    std::vector<std::vector<int> > adj_list;//[atom_index][adj_index]
     //[atom_index][min_distance between (atom with atom_index) and (atom with adj_index)]
    //std::vector<std::vector<vec3> > adj_dis_list;
    double mass;
    double Ek; //kinetic energy
    double temperature;
    
private:
    void randomv(double T_init, std::string v0_type);
    double v_fctr(double T_init);

};
#endif
