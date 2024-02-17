#ifndef __INTERSECTABLE_H__
#define __INTERSECTABLE_H__

#include "MathLib.h"
#include "Ray.h"
#include "Material.h"

using std::shared_ptr;
using std::make_shared;
using std::vector;

/******************************************************************
*                        Class Boundingbox                        *
******************************************************************/
class Boundingbox {
public:
    vec3 pmin, pmax;

    Boundingbox(): pmin(ZERO_VEC3), pmax(ZERO_VEC3) {}
    Boundingbox(vec3 pmin_, vec3 pmax_): 
    pmin(pmin_), pmax(pmax_) {}

    // Compute wheather a ray hitted the bounding box
    bool intersect(const Ray &ray) {
        float tmin = (pmin.x - ray.orig.x) / ray.dir.x;
        float tmax = (pmax.x - ray.orig.x) / ray.dir.x;

        if (tmin > tmax) std::swap(tmin, tmax);

        float tymin = (pmin.y - ray.orig.y) / ray.dir.y;
        float tymax = (pmax.y - ray.orig.y) / ray.dir.y;

        if (tymin > tymax) std::swap(tymin, tymax);

        if ((tmin > tymax) || (tymin > tmax)) return false;

        if (tymin > tmin) tmin = tymin;
        if (tymax < tmax) tmax = tymax;

        float tzmin = (pmin.z - ray.orig.z) / ray.dir.z;
        float tzmax = (pmax.z - ray.orig.z) / ray.dir.z;

        if (tzmin > tzmax) std::swap(tzmin, tzmax);

        if ((tmin > tzmax) || (tzmin > tmax)) return false;

        return true;
    }
};

/******************************************************************
*                        Class Interval                           *
******************************************************************/
class Intersectable {
public:
    Boundingbox bbox;
    virtual ~Intersectable() = default;
    // Compute wheather a ray hitted the object
    virtual bool intersect(const Ray &ray, Interval &rayt, InterRecord &rec) const = 0;
};

/******************************************************************
*                        Class Sphere                             *
******************************************************************/
class Sphere: public Intersectable {
public:
    double radius;  // radius of the sphere
    vec3 pos;  // position, emission, color of the sphere
    shared_ptr<Material> mat;

    // Sphere(vec3 p, double rad, shared_ptr<Material> mat): 
    // pos(p), radius(rad), mat(mat)) {}

    Sphere(vec3 p, double rad, shared_ptr<Material> mat) {
        pos = p;
        radius = rad;
        this->mat = mat;
        vec3 pmin = p - (XA + YA + ZA)*rad;
        vec3 pmax = p + (XA + YA + ZA)*rad;
        bbox = Boundingbox(pmin, pmax);
    }

    // Compute wheather a ray hitted the sphere
    bool intersect(const Ray &ray, Interval &rayt, InterRecord &rec) const override {
        // delta = (-b +- sqrt(b^2 - 4ac)) / 2
        // a = dir^2, b = 2*dir*op, c = op^2 - r^2
        // Assumed that ray.dir is normalized
        double t;
        vec3 op = pos - ray.orig;   // -op
        double b = ray.dir * op; // -b
        double det = b*b - (op*op - radius*radius);
        if (det < 0) return 0;  // onhit, return 0
        else det = std::sqrt(det);

        // return smaller root if it > 0, else return larger root if it > 0.
        // return 0 if both roots < 0.
        t = (t = b - det) > EPS ? t : ((t = b + det) > EPS ? t : 0);
        if (!rayt.surround(t)) return false;

        // Update InterRecord
        rec.t = t;
        rec.p = ray.at(t);
        rec.normal = (rec.p - pos) / radius;
        rec.inward = (ray.dir * rec.normal < 0);
        rec.mat = mat;
        return true;
    }
};

#endif