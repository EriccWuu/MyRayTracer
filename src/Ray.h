#ifndef __RAY_H__
#define __RAY_H__

#include "MathLib.h"

using std::shared_ptr;

class Material;

/******************************************************************
*                           Struct Ray                            *
******************************************************************/
struct Ray {
public:
    vec3 orig;  // Origin of the ray
    vec3 dir;   // Direction of the ray

    Ray() {}
    Ray(const vec3& origin, const vec3& direction) : orig(origin), dir(direction) {}

    vec3 origin() const  { return orig; }
    vec3 direction() const { return dir; }
    vec3 at(double t) const { return orig + t*dir; }    // Position at orig + t*dir
};

/******************************************************************
*                      Struct Inter_record                        *
******************************************************************/
struct InterRecord {
public:
    vec3 p;         // Ray position at orig + t*dir
    vec3 normal;    // Outward normal of intersectial face
    double t;       // Transmit time
    bool inward;    // Wheather the ray is going into objects
    shared_ptr<Material> mat;

    void set_normal(const Ray &ray, const vec3 &normal) {
        inward = ray.dir*normal < 0;
        this->normal = inward ? normal : -normal; 
    }
};

/******************************************************************
*                       Struct Interval                           *
******************************************************************/
struct Interval {
public:
    double min, max;

    Interval(): min(+INF), max(-INF) {} // Ddfault to be empty
    Interval(double min_, double max_): min(min_), max(max_) {}

    inline bool contain(double x) const { return min <= x && x <= max; }   // x in [min, max]
    inline bool surround(double x) const { return min < x && x < max; }    // x in (min, max)
    inline double clamp(double x) const { 
        // clamp x to [min, max]
        return (x < min) ? min : ((x > max) ? max : x);
    }

    static const Interval empty, universe;
};

const static Interval empty   (+INF, -INF);
const static Interval universe(-INF, +INF);

#endif