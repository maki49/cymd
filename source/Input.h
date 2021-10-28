#include <string>
class Input
{
public:
    Input();
    ~Input();
    void read_input_file(std::string& input_file);
    std::string geo_in_file;
    int natom;
    double rcut_potential;    //rcut for Leonard-Jones potential
    double extra_rcut_neighbor;   //extra distance for neighbor list 
    int max_neighbor;
    int whose_adj;  //which atom's adjacent atoms will be printed

    //parameters for L-J potential
    double epsilon;
    double sigma;

    double mass;
    int read_vel = 0;
    double init_temperature;  //(K)
    bool append = true;
};

