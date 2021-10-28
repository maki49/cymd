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
    }
    
}