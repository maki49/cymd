#ifndef _VEC3_
#define _VEC3_
#include <vector>
#include <cmath>

class vec3
{
private:

public:
    double x;
    double y;
    double z;
    double norm;

    vec3() {}
    vec3(double xin, double yin, double zin);
    ~vec3();
    vec3 operator+(vec3 a)
    {
        vec3 v(this->x + a.x,this->y + a.y,this->z + a.z);
        return v;
    }
    void operator+=(vec3 a)
    {
        this->x = this->x + a.x;
        this->y = this->y + a.y;
        this->z = this->z + a.z;
        this->norm = sqrt(this->x * this->x + this->y * this->y + this->z * this->z);
        return;
    }
    vec3 operator-(vec3 a)
    {
        vec3 v(this->x - a.x,this->y - a.y,this->z - a.z);
        return v;
    }
    vec3 operator-()
    {
        vec3 v(-this->x,-this->y,-this->z);
        return v;
    }
    void operator-=(vec3 a)
    {
        this->x = this->x - a.x;
        this->y = this->y - a.y;
        this->z = this->z - a.z;
        this->norm = sqrt(this->x * this->x + this->y * this->y + this->z * this->z);
        return;
    }
    vec3 operator*(double a)
    {
        vec3 v(this->x * a, this->y * a, this->z * a);
        return v;
    }
    void operator*=(double a)
    {
        this->x *= a;
        this->y *= a;
        this->z *= a;
        this->norm *= a;
        return;
    }
    vec3 operator/(double a)
    {
        vec3 v(this->x / a, this->y / a, this->z / a);
        return v;
    }
    vec3 vmodv(vec3 a)
    {
        vec3 v(fmod(this->x, a.x), fmod(this->y, a.y), fmod(this->z, a.z));
        return v;
    }
    vec3 vmodf(double a)
    {
        vec3 v(fmod(this->x, a), fmod(this->y, a), fmod(this->z, a));
        return v;
    }

};
#endif