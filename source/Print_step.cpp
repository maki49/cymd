#include "Print_step.h"
#include <fstream>
#include <iostream>
#include <iomanip>
Print_step::Print_step(){}
Print_step::Print_step(Geo& geo, LJ_pot& lj) :geo_step(geo), lj_step(lj){}
Print_step::~Print_step()
{
}

void Print_step::print_info(int istep)
{
    this->print_position(istep);
    this->print_velocity(istep);
    this->print_force(istep);
    this->print_log(istep);
}

void Print_step::print_position(int istep)
{
    //const double Bohr2A = 0.5291770;
    std::ofstream ofs;
    if(istep)
        ofs.open("position.txt", std::ios_base::app);
    else  
        ofs.open("position.txt");
    ofs << "STEP " << istep << ": ATOMIC_POSTION(Ang.) " << std::endl;
    ofs << std::setw(4) << "num" << std::setw(20) << "x(Ang.)" <<
        std::setw(20) << "y(Ang.)" << std::setw(20) << "z(Ang.)" << std::endl;
    for (int ia = 0; ia < this->geo_step.natom;++ia)
    {
        ofs << std::setprecision(geo_step.precision) << std::setw(4) << ia + 1
            << std::setw(20) << this->geo_step.r_t->at(ia).x
            << std::setw(20) << this->geo_step.r_t->at(ia).y
            << std::setw(20) << this->geo_step.r_t->at(ia).z << std::endl;
    }
    ofs << std::endl;
    ofs.close();
    return;
}

void Print_step::print_velocity(int istep)
{
    //const double Bohr2A = 0.5291770;
    std::ofstream ofs;
    if(istep)
        ofs.open("velocity.txt", std::ios_base::app);
    else  
        ofs.open("velocity.txt");
    ofs << "STEP " << istep << ": ATOMIC_VELOCITY(Ang./ps) " << std::endl;
    ofs << std::setw(4) << "num" << std::setw(20) << "vx(Ang./ps)" <<
        std::setw(20) << "vy(Ang./ps)" << std::setw(20) << "vz(Ang./ps)" << std::endl;
    for (int ia = 0; ia < this->geo_step.natom;++ia)
    {
        ofs << std::setprecision(geo_step.precision) << std::setw(4) << ia + 1
            << std::setw(20) << this->geo_step.atom_v[ia].x
            << std::setw(20) << this->geo_step.atom_v[ia].y
            << std::setw(20) << this->geo_step.atom_v[ia].z << std::endl;
    }
    ofs << std::endl;
    ofs.close();
    return;
}

void Print_step::print_force(int istep)
{
    //const double Bohr2A = 0.5291770;
    std::ofstream ofs;
    if(istep)
        ofs.open("force.txt", std::ios_base::app);
    else  
        ofs.open("force.txt");
    ofs << "STEP " << istep << ": ATOMIC_FORCE(eV/Ang.) " << std::endl;
    ofs << std::setw(4) << "num" << std::setw(20) << "Fx(eV/Ang.)" <<
        std::setw(20) << "Fy(eV/Ang.)" << std::setw(20) << "Fz(eV/Ang.)" << std::endl;
    for (int ia = 0; ia < this->geo_step.natom;++ia)
    {
        ofs << std::setprecision(geo_step.precision) << std::setw(4) << ia + 1
            << std::setw(20) << this->lj_step.atom_force[ia].x
            << std::setw(20) << this->lj_step.atom_force[ia].y
            << std::setw(20) << this->lj_step.atom_force[ia].z << std::endl;
    }
    ofs << std::endl;
    ofs.close();
    return;
}
void Print_step::print_log(int istep)
{
    //const double Bohr2A = 0.5291770;
    std::ofstream ofs;
    if(istep)
        ofs.open("run.log", std::ios_base::app);
    else  
        ofs.open("run.log");
    ofs << "STEP " << istep << std::endl;
    ofs << "Kinetic Energy (eV): " << std::setprecision(this->geo_step.precision)
        << this->geo_step.Ek<< std::endl;
    ofs << "Potential Energy (eV): " << std::setprecision(this->geo_step.precision)
        << this->lj_step.Ep << std::endl;
    ofs << "Total Energy (eV): " << std::setprecision(this->geo_step.precision)
        << this->geo_step.Ek + this->lj_step.Ep << std::endl;
    ofs << "Temperature (K): " << std::setprecision(this->geo_step.precision)
        << this->geo_step.temperature << std::endl;
    ofs << std::endl;
    ofs.close();
}