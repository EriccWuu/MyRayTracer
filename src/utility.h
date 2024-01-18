#ifndef UTILITY_H
#define UTILITY_H

#include "Matrix.h"

class Ray {
public:
    vec3 orig;
    vec3 dir;

    Ray() {}

    Ray(const vec3& origin, const vec3& direction) : orig(origin), dir(direction) {}

    vec3 origin() const  { return orig; }
    vec3 direction() const { return dir; }
    vec3 at(double t) const { return orig + t*dir; }
};

struct Sphere {
    double radius;
    vec3 pos, emis, color;

    Sphere(double rad, vec3 p, vec3 c, vec3 e={0,0,0}): 
    radius(rad), pos(p), emis(e), color(c) {}

    double intersect(const Ray &ray) const {
        // delta = (-b +- sqrt(b^2 - 4ac)) / 2
        // a = dir^2, b = 2*dir*op, c = op^2 - r^2
        vec3 op = pos - ray.orig;   // -op
        double t, eps = 1e-4;
        double b = ray.dir*op; // -b/2
        double det = b*b - op*op + radius*radius;
        if (det < 0) return 0;  // onhit, return 0
        else det = std::sqrt(det);

        // return smaller root if it > 0, else return larger root if it > 0.
        // return 0 if both roots < 0.
        return (t = b - det) > eps ? t : ((t = b + det) > eps ? t : 0);
    }
};

#endif