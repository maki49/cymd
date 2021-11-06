#include "Input.h"
#include <iostream>
#include <fstream>
#include <string>

Input::Input(){}

Input::~Input() {}

void Input::read_input_file(std::string& input_file)
{
    std::ifstream ifs;
    ifs.open(input_file);
    if (!ifs)
    {
        std::cout << "Error in reading input file !" << std::endl;
        return;
	}
    ifs.clear();
    ifs.seekg(0);   //Sets the position of the next character
    ifs.rdstate();
    std::string tmp;
    while (ifs.good())
    {
        ifs >> tmp;
        if (tmp == "natom") ifs >> this->natom;
        if (tmp == "geo_in_file") ifs >> this->geo_in_file;
        if (tmp == "rcut_potential") ifs >> this->rcut_potential;
        if (tmp == "extra_rcut_neighbor") ifs >> this->extra_rcut_neighbor;
        if (tmp == "max_neighbor") ifs >> this->max_neighbor;
        if (tmp == "whose_adj") ifs >> this->whose_adj;
        if (tmp == "epsilon") ifs >> this->epsilon;
        if (tmp == "sigma") ifs >> this->sigma;
        if (tmp == "mass") ifs >> this->mass;
        if (tmp == "read_vel") ifs >> this->read_vel;
        if (tmp == "init_temperature") ifs >> this->init_temperature;
        if (tmp == "append") ifs >> this->append;
        if (tmp == "nstep") ifs >> this->nstep;
        if (tmp == "dt") ifs >> this->dt;
        if (tmp == "steps_per_print") ifs >> this->steps_per_print;
        if (tmp == "steps_per_search") ifs >> this->steps_per_search;
        if (tmp == "ensemble") ifs >> this->ensemble;
        if (tmp == "verlet_method") ifs >> this->verlet_method;
    }
    
}