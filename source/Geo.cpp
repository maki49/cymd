#include <cstdio>
#include "Geo.h"
#include <fstream>
#include <cstring>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <assert.h>
#include "Input.h"
#include <random>
#include <time.h>

Geo::Geo(){}
Geo::Geo(double m, double T0) : mass(m), temperature(T0) {}
Geo::~Geo() {}
void Geo::read_geo(std::string& geo_in_file, int read_vel)
{
    std::ifstream ifs;
    ifs.open(geo_in_file);
    if (!ifs)
    {
        std::cout << "Error in reading geo file !" << std::endl;
        return;
	}
    ifs.clear();
    //ifs.seekg(0);   //Sets the position of the next character
    //ifs.rdstate();
    std::string type_read;
    std::string tmp_str;
    double x_read;
    double y_read;
    double z_read;
    std::string x_str;

    ifs >> tmp_str;
    assert(tmp_str == "%CELL_PARAMETER");
    //read cell parameter
    for (int i = 0;i < 3;++i)
    {
        ifs >> x_read >> y_read >> z_read;
        vec3 v(x_read, y_read, z_read);
        this->lattice_vec[i] = v;
    }
    this->R = this->lattice_vec[0] + lattice_vec[1] + lattice_vec[2];

    ifs >> tmp_str;
    while (ifs.good() && tmp_str != "%ATOMIC_POSTION")
    {
        ifs >> tmp_str;
    }
    //set precision according to the first data
    if (ifs.good() && this->precision == 0)
    {
        ifs >> type_read >> x_str >> y_read >> z_read;
        this->precision = x_str.length() - 1;
        /*
        //for fixed precision
        for (int i = 0;i < x_str.length();++i)
        {
            if (x_str[i] == '.')
            {
                this->precision = x_str.length() - i - 1;
                break;
            }
        }
        */
        std::stringstream ss;
        ss << x_str;
        ss >> x_read;
        vec3 vec_read(x_read, y_read, z_read);
        this->atom_types.push_back(type_read);
        this->atom_coords.push_back(vec_read);
    }

    while (ifs.good())
    {
        ifs >> type_read;
        //find precision
        if (type_read == "%ATOMIC_VELOCITY") break;
        ifs >> x_read >> y_read >> z_read;
        this->atom_types.push_back(type_read);
        vec3 vec_read(x_read, y_read, z_read);
        this->atom_coords.push_back(vec_read);
    }
    this->natom = this->atom_coords.size();
    this->adj_list.resize(this->natom);
    this->adj_dis_list.resize(this->natom);
    
    // when ATOMIC_VELOCITY is not needed
    

    if (!read_vel)  //generate velo radomly
    {
        ifs.close();
        this->randomv(this->temperature);
        std::cout << "Generated velocity randomly." << std::endl;
        return;
    }
    else
    {
        if (type_read != "%ATOMIC_VELOCITY")
        {
            std::cout << "Please provide '%ATOMIC_VELOCITY' in `geo_file` when `read_val` is 1." << std::endl;
            exit(0);
        }
    }
    //continue reading velocity
    double sum_v2 = 0;
    for (int i = 0;i < this->natom;++i)
    {
        ifs >> type_read;
        if (type_read == "")break;
        //find precision
        ifs >> x_read >> y_read >> z_read;
        vec3 vec_read(x_read, y_read, z_read);
        this->atom_v.push_back(vec_read);

        sum_v2 += vec_read.norm * vec_read.norm;
    }
    assert(this->natom == this->atom_v.size());
    std::cout << "Read atom velocity from geo file." << std::endl;
    this->cal_EkT(sum_v2);
    ifs.close();
    return;
    // remain: same precision
}

void Geo::search_adj(double rcut, int max_neighbor)
{
    for (int i = 0;i < this->natom;++i)
    {
        this->adj_dis_list[i].clear();
        this->adj_list[i].clear();
    }
    for (int i = 0;i < this->natom;++i)
    {
        for (int j = i;j < this->natom;++j)
        {
            vec3 nearest_dr = (this->atom_coords[i] - this->atom_coords[j]).vmodv(this->R);
            double min_distance = nearest_dr.norm;
            //if (i == 11) std::cout << j + 1 << " " << 0 << " " <<0<< " " <<0<<" "<<min_distance<<std::endl;
            for (int fx = -1;fx <= 1;++fx)
            {
                for (int fy = -1;fy <= 1;++fy)
                {
                    for (int fz = -1;fz <= 1;++fz)
                    {
                        vec3 dR(R.x * (double)fx, R.y * (double)fy, R.z * (double)fz);
                        vec3 dr = (this->atom_coords[i] - this->atom_coords[j] + dR).vmodv(R);  //caution: ri-rj
                        if (dr.norm < min_distance)
                        {
                            min_distance = dr.norm;
                            nearest_dr = dr;
                            //if (i == 11) std::cout << j + 1 << " " << fx << " " <<fy<< " " <<fz<<" "<<min_distance<<std::endl;
                        }
                    }
                }
            }
            if (min_distance <= rcut && this->adj_list[i].size() <= max_neighbor)
            {
                if (j == i) continue;  //exclude itself
                //including itself
                this->adj_list[i].push_back(j);
                this->adj_dis_list[i].push_back(nearest_dr);
                if (i != j)
                {
                    this->adj_list[j].push_back(i);
                    this->adj_dis_list[j].push_back(-nearest_dr);   //caution: rji = - rij
                }
                //assert (nearest_dr.norm == min_distance);
                assert (abs(nearest_dr.norm - min_distance)<1e-10);
            }
        }
    }
    return;
}

void Geo::print_adj_list(int ia)
{
    const double Bohr2A = 0.5291770;
    std::ofstream ofs;
    ofs.open("geo.out");
    ofs << "the neighbor list of atom " << ia << std::endl;
    ofs << std::setw(8) << "Num " << std::setw(15) << "x(Ang)"<< std::setw(15) << "y(Ang)" << std::setw(15) << "z(Ang)" << std::endl;
    for (int j = 0;j < this->adj_list[ia - 1].size();++j)
    {
        if (this->adj_list[ia - 1][j] + 1 == ia) continue;
        ofs << std::setw(8) << this->adj_list[ia - 1][j] + 1 << " "
            //<< std::setw(15) << std::setprecision(this->precision) << this->adj_dis_list[ia - 1][j]
            //<< std::setw(15) << std::setprecision(this->precision) << this->adj_dis_list[ia - 1][j]/Bohr2A << std::endl;
            << std::setw(15) << std::setprecision(this->precision) << this->atom_coords[this->adj_list[ia-1][j]].x <<" "
            << std::setw(15) << std::setprecision(this->precision) << this->atom_coords[this->adj_list[ia-1][j]].y <<" "
            << std::setw(15) << std::setprecision(this->precision) << this->atom_coords[this->adj_list[ia-1][j]].z <<std::endl;
    }
    ofs.close();
}

void Geo::randomv(double T_init)
{
    //step1: generate randomly and sum
    srand((unsigned)time(NULL));    //random seed
    
    vec3 sum_v = vec3(0, 0, 0);
    for (int i = 0;i < this->natom;++i)
    {
        vec3 v = vec3(rand() / double(RAND_MAX) - 0.5,
            rand() / double(RAND_MAX) - 0.5, rand() / double(RAND_MAX) - 0.5);
        this->atom_v.push_back(v);
        //sum_mv+=v*atom_mass[atom_type[ia]];
        //sum_m+=atom_mass[atom_type[ia]];
        sum_v+=v;   // for all the atoms are the same type
    }

    //step2: v -> (v-vc)
    //vc = sum_mv/sum_m;     // velocity of mass centor
    vec3 vc = sum_v / double(this->natom);      // velocity of mass centor, for same atoms
    for (int i = 0;i < this->natom;++i)
    {
        this->atom_v[i] -= vc;
    }

    //3. factorize
    double sum_v2 = 0;
    double vf = v_fctr(T_init);
    for (int i = 0;i < this->natom;++i)
    {
        this->atom_v[i] *= vf;
        sum_v2 += this->atom_v[i].norm * this->atom_v[i].norm;
    }
    assert(this->atom_v.size() == natom);
    
    //4. cal Ek and T
    this->cal_EkT(sum_v2);
    return;
}

double Geo::v_fctr(double T_init)
{
    //for all atoms are in the same type
    double sum_v2 =0.0;
    for (auto v : this->atom_v)
    {
        sum_v2 += v.norm * v.norm;
    }
    return sqrt(3 * double(this->natom) * kb * T_init / (mass / NA * 1e-3) / (sum_v2 * 1e4));
}

void Geo::cal_EkT(double sum_v2)
{
    //cal temperature by 
    //1e4: (Ang/ps)^2 -> (m/s)^2
    this->temperature = (mass / NA *1e-3) * (sum_v2 * 1e4) / 3 / double(this->natom) / kb;
    this->Ek = 1.5 * double(this->natom) * kb * this->temperature / e *1e-4 ; //e: J -> eV; -4=-23+19
    this->Ek = 0.5 * (mass / NA * 1e-3) * (sum_v2 * 1e4) / e * 1e-4;  //1e(23-19)
    return;
}

void Geo::search_adj_faster(double rcut, int max_neighbor)
{
    for (int i = 0;i < this->natom;++i)
    {
        this->adj_dis_list[i].clear();
        this->adj_list[i].clear();
    }
    for (int i = 0;i < this->natom;++i)
    {
        for (int j = i + 1;j < this->natom;++j)
        {
            double L = this->lattice_vec[0].norm;
            vec3 nearest_dr = this->shortest(this->atom_coords[i] - this->atom_coords[j]);
            double min_distance = nearest_dr.norm;
            //judge if is adjacent atom
            if (min_distance <= rcut && this->adj_list[i].size() <= max_neighbor)
            {
                if (j == i) continue;  //exclude itself
                //including itself
                this->adj_list[i].push_back(j);
                this->adj_dis_list[i].push_back(nearest_dr);
                if (i != j)
                {
                    this->adj_list[j].push_back(i);
                    this->adj_dis_list[j].push_back(-nearest_dr);   //caution: rji = - rij
                }
                //assert (nearest_dr.norm == min_distance);
                assert (abs(nearest_dr.norm - min_distance)<1e-10);
            }
        }
    }
    return;
}
void Geo::update_dis_list(void)
{
    for (int i = 0;i < this->natom;++i)
    {
        int index_j = 0;
        for (auto j : this->adj_list[i])
        {
            //means all elements after j in adj_list[i] >= i, 
            //and they will be calculated at adj_list[j] and so on
            vec3 dr = this->shortest(this->atom_coords[i] - this->atom_coords[j]);
            this->adj_dis_list[i][index_j] = dr;
            index_j += 1;
        }
    }
    return;
}

vec3 Geo::res_in_box(vec3 r)
{
    r = r.vmodv(this->R);
    if (r.x < 0)r.x += this->R.x;
    if (r.y < 0)r.y += this->R.y;
    if (r.z < 0)r.z += this->R.z;
    return r;
}
//res_in_box(dr+R/2)-R/2 is same as "shortest(dr, R)"
vec3 Geo::shortest(vec3 r)
{
    return this->res_in_box(r + this->R / 2) - this->R / 2;
}