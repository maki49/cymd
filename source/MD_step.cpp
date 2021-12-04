#include "MD_step.h"
#include "iostream"
#include <time.h>
#include <cmath>
#include <iomanip>
#include <assert.h>

MD_step::MD_step() {}
MD_step::MD_step(Input& input) :
    dt(input.dt),
    nstep(input.nstep),
    steps_per_print(input.steps_per_print),
    steps_per_search(input.steps_per_search),
    rcut_neighbor(input.rcut_potential + input.extra_rcut_neighbor),
    max_neighbor(input.max_neighbor),
    verlet_method(input.verlet_method),
    append(input.append),
    ensemble(input.ensemble),
    cal_msd(input.cal_msd),
    msd_print_interval(input.msd_print_interval),
    cal_rdf(input.cal_rdf),
    rdf_rcut(input.rdf_rcut),
    rdf_dr(input.rdf_dr),
    rdf_start_step(input.rdf_start_step),
    rdf_end_step(input.rdf_end_step),
    rdf_interval(input.rdf_interval)
{
    if (ensemble == "NVT")
    {
        this->thermo_temperature = input.thermo_temperature;
        this->nraise = input.nraise;
    }
};
MD_step::~MD_step() {}

void MD_step::allocate(int natom)
{
    switch (this->verlet_method)
    {
    case 0:
        this->r1.resize(natom);
        std::cout << "Algorithm: Verlet (method 1)" << std::endl;
        break;
    case 1:
        this->r1.resize(natom);
        this->r2.resize(natom);
        std::cout << "Algorithm: Verlet (method 2)" << std::endl;
        break;
    case 2:
        std::cout << "Algorithm: Velocity Verlet" << std::endl;
        break;
    }

    if (this->cal_msd)
    {
        vec3 v0(0,0,0);
        this->atom_msd.resize(natom, v0);

        std::ofstream ofs;
        ofs.open("msd.txt");
        ofs<<std::setw(4)<<"t/ps"<<std::setw(20)<<"msd (Ang.^2)"<< std::endl;
        ofs.close();
    }

    if (this->cal_rdf)
    {
        assert(this->rdf_dr > 1e-5);
        this->nshell_rdf.resize((int)(this->rdf_rcut / this->rdf_dr) + 1, 0.0);
    }
        

    return;
}

void MD_step::verlet_1(int istep, Geo& geo_step, LJ_pot& lj_step)
{
    double dt2 = dt * dt;
    //double m = geo_step.mass * 1e-3 / geo_step.e / geo_step.NA; // g/mol -> eV*ps^2/Ang.^2
    if (istep == 0)
    {
        for (int ia = 0;ia < geo_step.natom;++ia)
        {
            this->r_tpdt->at(ia) =
                geo_step.res_in_box(geo_step.r_t->at(ia)
                + geo_step.atom_v[ia] * dt
                + lj_step.atom_force[ia] * 0.5 / m * dt2, geo_step.R);
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
            vec3 v = geo_step.shortest(
                this->r_tpdt->at(ia) - this->r_tmdt->at(ia), geo_step.R) / 2 / dt;
            sum_v2 += v.norm * v.norm;  //(Ang/ps)^2
            //nvc += v;
            geo_step.atom_v[ia] = v;
            r_tpdt->at(ia) = geo_step.res_in_box(r_tpdt->at(ia), geo_step.R);
        }
        //std::cout <<nvc.norm << std::endl;
        geo_step.cal_EkT(sum_v2);
    }
    return;
}

void MD_step::verlet_0(int istep, Geo& geo_step, LJ_pot& lj_step)
{
    double dt2 = dt * dt;
    if (istep == 0)
    {
        for (int ia = 0;ia < geo_step.natom;++ia)
        {
            this->r_tmdt->at(ia) =//actually tpdt
                geo_step.res_in_box(geo_step.r_t->at(ia)
                + geo_step.atom_v[ia] * dt
                    + lj_step.atom_force[ia] * 0.5 / m * dt2,
                    geo_step.R);
        }
    }
    else
    {
        double sum_v2 = 0.0;
        for (int ia = 0;ia < geo_step.natom;++ia)
        {
            vec3 v = geo_step.shortest(
                geo_step.r_t->at(ia) - this->r_tmdt->at(ia), geo_step.R) / dt
                + lj_step.atom_force[ia] / m / 2 * dt;
            sum_v2 += v.norm * v.norm;  //(Ang/ps)^2
            geo_step.atom_v[ia] = v;
            this->r_tmdt->at(ia) = geo_step.res_in_box(//actually tpdt
                geo_step.r_t->at(ia) * 2 - this->r_tmdt->at(ia)
                + lj_step.atom_force[ia] / m * dt2, geo_step.R);
        }
        geo_step.cal_EkT(sum_v2);
    }
    return;
}

//week 8
void MD_step::velocity_verlet_before(Geo& geo_step, LJ_pot& lj_step)
{
    double dt2 = dt * dt;
    for (int ia = 0;ia < geo_step.natom;++ia)
    {
        vec3 dr = geo_step.atom_v[ia] * dt
            + lj_step.atom_force[ia] / 2 / m * dt2;
        geo_step.r_t->at(ia) += dr;
        if(this->cal_msd)
            this->atom_msd.at(ia)+=geo_step.shortest(dr, geo_step.R);
        geo_step.atom_v[ia] += lj_step.atom_force[ia] / 2 / m * dt;
    }
    return;
}
void MD_step::velocity_verlet_after(Geo& geo_step, LJ_pot& lj_step)
{
    double sum_v2 = 0.0;
    for (int ia = 0;ia < geo_step.natom;++ia)
    {
        geo_step.atom_v[ia] += lj_step.atom_force[ia] / 2 / m * dt;
        sum_v2 += geo_step.atom_v[ia].norm * geo_step.atom_v[ia].norm;
    }
    geo_step.cal_EkT(sum_v2);
    return;
}

void MD_step::Anderson(Geo& geo_step, double sgm, 
    std::default_random_engine& generator)
{
    auto gaussrand = [&]() ->double
    {   //std-lib
        std::normal_distribution<double> distribution(0, sgm);
        return distribution(generator);
    };

    auto gaussrand_bm = [&]() ->double
    {   //Box-Muller Algorithm
        double u1 = rand() / double(RAND_MAX);
        double u2 = rand() / double(RAND_MAX);
        return sqrt(-2 * log(u1)) * sin(2 * M_PI * u2) * sgm;
    };

    auto gaussrand_1 = [&]() ->double
    {
        double S=0.0;
        double v1, v2, u1, u2;
        while (S>=1||S==0.0){
            u1 = rand() / double(RAND_MAX);
            u2 = rand() / double(RAND_MAX);
            v1=2*u1-1;
            v2=2*u2-1;
            S=v1*v1+v2*v2;
        }
        double X=v1*sqrt(-2*log(S)/S);
        return X*sgm;
    };

    for (int i = 0;i < geo_step.natom;++i)
    {
        if ((double(rand()) / double(RAND_MAX)) < 1.0 / nraise)
        {
            vec3 v(gaussrand(), gaussrand(), gaussrand()); 
            //vec3 v(gaussrand_bm(), gaussrand_bm(), gaussrand_bm()); 
            //vec3 v(gaussrand_1(), gaussrand_1(), gaussrand_1()); 
            geo_step.atom_v[i] = v;
        }
    }
    return;
}

void MD_step::main_step(Geo& geo_step, LJ_pot& lj_step)
{
    this->allocate(geo_step.natom); //allocate memory 
    //mass: g/mol -> eV*ps^2/Ang.^2
    this->m = geo_step.mass * 1e-3 / geo_step.e / geo_step.NA;
    //set sigma for NVE
    double sgm_ads = sqrt(geo_step.kb * this->thermo_temperature
        / (geo_step.mass * 1e-3 / geo_step.NA))*1e-2; 

    //random seed
    //( Caution: FAKE random number, set seed INITIALLY !! )
    srand((unsigned)time(NULL));
    std::default_random_engine generator(time(NULL));

    std::cout << "===== MD start ======" << std::endl;
    clock_t start, end;
    start = clock();
    
    for (int istep = 0; istep <= this->nstep;++istep)
    {
        //0. for velocity verlet: update r_t and v_t before f_t
        if (this->verlet_method == 2 && istep > 0)
            this->velocity_verlet_before(geo_step, lj_step);

        //1. update adj_list of current step (with r_t)
        if (istep % this->steps_per_search == 0)
        {
            geo_step.search_adj_faster(this->rcut_neighbor, max_neighbor);
            //geo_step.print_adj_list(12);
        }
        //else//only update adj_dis_list
            //geo_step.update_dis_list();
        
        //2. cal force(f_t) based on geo of current step
        lj_step.cal_EpF(geo_step);

        //3. verlet: calculate r_t+dt, v_t and Ek;
        //velocity verlet: update v_t after f_t, and calculate Ek
        switch (this->verlet_method)
        {
        case 0:
            this->verlet_0(istep, geo_step, lj_step);
            break;
        case 1:
            this->verlet_1(istep, geo_step, lj_step);
            break;
        case 2:
            if (istep > 0)
                this->velocity_verlet_after(geo_step, lj_step);
            break;
        }
        //add velo-verlet here
        
        //4. print info(r,v,f) of current step 
        Print_step ps(geo_step, lj_step);
        if (istep % this->steps_per_print == 0)
            ps.print_info(istep, this->append);

        //5. update geo(r_t+dt) for next step 
        //(exchange memory by pointers)
        //no need in velocity verlet
        std::vector<vec3>* tmp;
        tmp = this->r_tmdt;
        switch (this->verlet_method)
        {
        case 0:
            this->r_tmdt = geo_step.r_t;
            geo_step.r_t = tmp;
            break;
        case 1:
            this->r_tmdt = geo_step.r_t;
            geo_step.r_t = this->r_tpdt;
            this->r_tpdt = tmp;
            break;
        }

        // 6. collision, if NVT
        if (this->ensemble == "NVT")
            this->Anderson(geo_step, sgm_ads, generator);

        //7. msd(if needed)
        if (this->cal_msd && istep%this->msd_print_interval==0)
            this->print_msd(istep, geo_step.natom, geo_step.precision);

        //8. rdf(if needed)
        if (this->cal_rdf && istep >= this->rdf_start_step &&
            istep <= this->rdf_end_step && istep % rdf_interval == 0)
            this->rdf(geo_step);
    }
    
    end = clock();
    std::cout << "=====================" << std::endl;
    std::cout << "MD finished in  " << double(end - start) / CLOCKS_PER_SEC
        << " seconds for " << this->nstep << " steps." << std::endl;
    
    if (this->cal_rdf)
        this->print_rdf((double)geo_step.natom
            / (geo_step.R.x * geo_step.R.y * geo_step.R.z));

    return;
}


void MD_step::print_msd(int istep, int natom, int precision)
{
    std::ofstream ofs;
    ofs.open("msd.txt", std::ios_base::app);
    double sum_msd = 0.0;
    for (int ia = 0; ia < natom;++ia)
    {
        sum_msd += this->atom_msd[ia].norm * this->atom_msd[ia].norm;
    }
    ofs << std::setw(4) << istep * this->dt << std::setw(20)
        << std::setprecision(precision) << sum_msd / natom << std::endl;
    ofs.close();
    return;
}

void MD_step::rdf(Geo& geo_step)
{
    //average of each atom
    double invnatom = 1.0 / double(geo_step.natom);
    //if rdf_rcut < rcut_adj
    for (int ia = 0;ia < geo_step.natom;++ia)
    {
        for (auto ja : geo_step.adj_list[ia])
        {
            vec3 rij = Geo::shortest(geo_step.r_t->at(ia)
                - geo_step.r_t->at(ja), geo_step.R);
            if (rij.norm < this->rdf_rcut)
                this->nshell_rdf.at((int)(rij.norm / this->rdf_dr))+=invnatom;
        }
    }
    this->rdf_ncal += 1;
    return;
}

void MD_step::print_rdf(double rho)
{
    std::cout << rho << std::endl;
    std::ofstream ofs;
    ofs.open("rdf.txt");
    for (int ir = 0;ir < this->nshell_rdf.size()-1;++ir)
    {
        double r = (double)ir * this->rdf_dr + 0.5 * this->rdf_dr;
        double dv = 4 * M_PI * r * r * this->rdf_dr;
        ofs << std::setw(4) << std::setprecision(2) << r
            << std::setw(20) << std::setprecision(6)
            << (double)nshell_rdf.at(ir) / dv / rho
            / (double)this->rdf_ncal << std::endl;
    }
    ofs.close();
    return;
}
