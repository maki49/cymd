#include "vec3.h"
#include <cmath>

vec3::vec3(double xin, double yin, double zin) : x(xin), y(yin), z(zin)
{
    this->norm = sqrt(this->x * this->x + this->y * this->y + this->z * this->z);
}

vec3::~vec3(){}
