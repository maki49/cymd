#include "MD_step.h"
#include "iostream"
MD_step::MD_step() {}
MD_step::MD_step(Input& input) :
    dt(input.dt),
    nstep(input.nstep),
    steps_per_print(input.steps_per_print),
    steps_per_search(input.steps_per_search),
    rcut_neighbor(input.rcut_potential + input.extra_rcut_neighbor),
    verlet_method(input.verlet_method)
{};
MD_step::~MD_step() {}

void MD_step::allocate(int natom)
{
    switch (this->verlet_method)
    {
    case 1:
        this->r1.resize(natom);
        break;
    case 2:
        this->r1.resize(natom);
        this->r2.resize(natom);
        break;
    }
    return;
}

void MD_step::verlet_2(int istep, Geo& geo_step, LJ_pot& lj_step)
{
    double dt2 = dt * dt;
    double m = geo_step.mass * 1e-3 / geo_step.e / geo_step.NA; // g/mol -> eV*ps^2/Ang.^2
    if (istep == 0)
    {
        for (int ia = 0;ia < geo_step.natom;++ia)
        {
            this->r_tpdt->at(ia) =
                geo_step.res_in_box(geo_step.r_t->at(ia)
                + geo_step.atom_v[ia] * dt
                + lj_step.atom_force[ia] * 0.5 / m * dt2);
        }
    }
    else
    {
        double sum_v2 = 0.0;
        //vec3 nvc(0, 0, 0);
        for (int ia = 0;ia < geo_step.natom;++ia)
        {
            this->r_tpdt->at(ia) = geo_step.r_t->at(ia) * 2
                - this->r_tmdt->at(ia) + lj_step.atom_force[ia] / m * dt2;
            vec3 v = geo_step.shortest(this->r_tpdt->at(ia) - this->r_tmdt->at(ia)) / 2 / dt;
            sum_v2 += v.norm * v.norm;  //(Ang/ps)^2
            //nvc += v;
            geo_step.atom_v[ia] = v;
            r_tpdt->at(ia) = geo_step.res_in_box(r_tpdt->at(ia));
        }
        //std::cout <<nvc.norm << std::endl;
        geo_step.cal_EkT(sum_v2);
    }
    return;
}

void MD_step::verlet_1(int istep, Geo& geo_step, LJ_pot& lj_step)
{
    double dt2 = dt * dt;
    double m = geo_step.mass * 1e-3 / geo_step.e / geo_step.NA; // g/mol -> eV*ps^2/Ang.^2
    if (istep == 0)
    {
        for (int ia = 0;ia < geo_step.natom;++ia)
        {
            this->r_tmdt->at(ia) =//actually tpdt
                geo_step.res_in_box(geo_step.r_t->at(ia)
                + geo_step.atom_v[ia] * dt
                + lj_step.atom_force[ia] * 0.5 / m * dt2);
        }
    }
    else
    {
        double sum_v2 = 0.0;
        for (int ia = 0;ia < geo_step.natom;++ia)
        {
            vec3 v = geo_step.shortest(geo_step.r_t->at(ia) - this->r_tmdt->at(ia)) / dt
                + lj_step.atom_force[ia] / m / 2 * dt;
            sum_v2 += v.norm * v.norm;  //(Ang/ps)^2
            geo_step.atom_v[ia] = v;
            this->r_tmdt->at(ia) = geo_step.res_in_box(//actually tpdt
                geo_step.r_t->at(ia) * 2 - this->r_tmdt->at(ia)
                + lj_step.atom_force[ia] / m * dt2);
        }
        geo_step.cal_EkT(sum_v2);
    }
    return;
}

//week 8
void MD_step::velocity_verlet(int istep, Geo& geo_step, LJ_pot& lj_step)
{
    
}

void MD_step::main_step(Geo& geo_step, LJ_pot& lj_step)
{
    this->allocate(geo_step.natom); //allocate memory 
    for (int istep = 0; istep <= this->nstep;++istep)
    {
        //1. update adj_list of current step (with r_t)
        if (istep % this->steps_per_search == 0)
            geo_step.search_adj_faster(this->rcut_neighbor, max_neighbor);
        else//only update adj_dis_list
            geo_step.update_dis_list();
        
        //2. cal force(f_t) based on geo of current step
        lj_step.cal_EpF(geo_step);

        //3. verlet: calculate r_t+dt, v_t and Ek
        switch (this->verlet_method)
        {
        case 1:
            this->verlet_1(istep, geo_step, lj_step);
            break;
        case 2:
            this->verlet_2(istep, geo_step, lj_step);
            break;
        }
        //add velo-verlet here
        
        //4. print info(r,v,f) of current step 
        Print_step ps(geo_step, lj_step);
        if (istep % this->steps_per_print == 0)
            ps.print_info(istep, 1);

        //5. update geo(r_t+dt) for next step (exchange memory by pointers)
        std::vector<vec3>* tmp;
        tmp = this->r_tmdt;
        switch (this->verlet_method)
        {
        case 1:
            this->r_tmdt = geo_step.r_t;
            geo_step.r_t = tmp;
            break;
        case 2:
            this->r_tmdt = geo_step.r_t;
            geo_step.r_t = this->r_tpdt;
            this->r_tpdt = tmp;
            break;
        }
        
    }
    return;
}
