#ifndef UTILITY_H
#define UTILITY_H

#include "Matrix.h"

class ray {
public:
    ray() {}

    ray(const vec3& origin, const vec3& direction) : orig(origin), dir(direction) {}

    vec3 origin() const  { return orig; }
    vec3 direction() const { return dir; }

    vec3 at(double t) const {
        return orig + t*dir;
    }
    vec3 tmp(vec3 v) ;

private:
    vec3 orig;
    vec3 dir;
};

#endif