#ifndef UTILITY_H
#define UTILITY_H

#include "MathLib.h"
#include <memory>
#include <vector>


// Usings
using std::shared_ptr;
using std::make_shared;
using std::vector;

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
struct Inter_record {
public:
    vec3 p;         // Ray position at orig + t*dir
    vec3 normal;    // Outward normal of intersectial face
    double t;       // Transmit time
    bool inward;    // Wheather the ray is going into objects

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
        return (x < min) ? min : ((min > x) ? max : x);
    }

    static const Interval empty, universe;
};

const static Interval empty   (+INF, -INF);
const static Interval universe(-INF, +INF);

/******************************************************************
*                        Class Interval                           *
******************************************************************/
class Intersectable {
public:
    virtual ~Intersectable() = default;
    // Compute wheather a ray hitted the object
    virtual bool intersect(const Ray &ray, Interval &rayt, Inter_record &rec) const = 0;
};

typedef std::vector<std::shared_ptr<Intersectable>> Interlist;

/******************************************************************
*                        Class Sphere                             *
******************************************************************/
class Sphere: public Intersectable {
public:
    double radius;  // radius of the sphere
    vec3 pos, emis, color;  // position, emission, color of the sphere

    Sphere(vec3 p, double rad, vec3 c={1, 1, 1}, vec3 e={0,0,0}): 
    pos(p), radius(rad), color(c), emis(e) {}

    // Compute wheather a ray hitted the sphere
    bool intersect(const Ray &ray, Interval &rayt, Inter_record &rec) const override {
        // delta = (-b +- sqrt(b^2 - 4ac)) / 2
        // a = dir^2, b = 2*dir*op, c = op^2 - r^2
        double t, eps = 1e-4;
        vec3 op = pos - ray.orig;   // -op
        double a_recip = 1 / ray.dir.norm2();
        double b = ray.dir * op * a_recip; // -b/2a
        double det = b*b - (op*op - radius*radius)*a_recip;
        if (det < 0) return 0;  // onhit, return 0
        else det = std::sqrt(det);

        // return smaller root if it > 0, else return larger root if it > 0.
        // return 0 if both roots < 0.
        t = (t = b - det) > eps ? t : ((t = b + det) > eps ? t : 0);
        if (!rayt.surround(t)) return false;
        rec.t = t;
        rec.p = ray.at(t);
        rec.set_normal(ray, (rec.p - pos) / radius);
        return true;
    }
};

#endif