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
class BoundingBox {
public:
    vec3 pmin, pmax;

    BoundingBox(): pmin(ZERO_VEC3), pmax(ZERO_VEC3) {}
    BoundingBox(const vec3 &pmin_, const vec3 &pmax_): 
    pmin(pmin_), pmax(pmax_) {}
    BoundingBox(const BoundingBox &a, const BoundingBox &b) {
        float ax = std::min(a.pmin.x, b.pmin.x);
        float ay = std::min(a.pmin.y, b.pmin.y);
        float az = std::min(a.pmin.z, b.pmin.z);
        float bx = std::max(a.pmax.x, b.pmax.x);
        float by = std::max(a.pmax.y, b.pmax.y);
        float bz = std::max(a.pmax.z, b.pmax.z);
        pmin = vec3(ax, ay, az);
        pmax = vec3(bx, by, bz);
    }
    // Compute wheather a ray hitted the bounding box
    bool intersect(const Ray &ray) const {
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
*                        Class Intersectable                      *
******************************************************************/
class Intersectable {
public:
    BoundingBox bbox;
    virtual ~Intersectable() = default;
    // Compute wheather a ray hitted the object
    virtual bool intersect(const Ray &ray, Interval rayt, InterRecord &rec) const = 0;
    virtual BoundingBox boundingBox() const = 0;
};

static bool compare(const shared_ptr<Intersectable> &a, const shared_ptr<Intersectable> &b, int axis = 0) {
    return a->boundingBox().pmin[axis] < b->boundingBox().pmin[axis];
}

/******************************************************************
*                        Class Sphere                             *
******************************************************************/
class Sphere: public Intersectable {
public:
    double radius;  // radius of the sphere
    vec3 pos;  // position, emission, color of the sphere
    shared_ptr<Material> mat;

    Sphere(vec3 p, double rad, shared_ptr<Material> mat) {
        pos = p;
        radius = rad;
        this->mat = mat;
        bbox = boundingBox();
    }

    BoundingBox boundingBox() const override {
        vec3 pmin = pos - (XA + YA + ZA)*radius;
        vec3 pmax = pos + (XA + YA + ZA)*radius;
        return BoundingBox(pmin, pmax);
    }

    static vec2 uv(const vec3 &p) {
        double theta = acos(-p.y);
        double phi = atan2(-p.z, p.x) + PI;
        return vec2(phi / 2 / PI, theta / PI);
    }

    // Compute wheather a ray hitted the sphere
    bool intersect(const Ray &ray, Interval rayt, InterRecord &rec) const override {
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
        rec.uv = uv(rec.normal);
        rec.inward = (ray.dir * rec.normal < 0);
        rec.mat = mat;

        return true;
    }

};

#endif