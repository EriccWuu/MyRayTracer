#ifndef __SPHERE_H__
#define __SPHERE_H__
// #include "MathLib.h"
#include "Intersectable.h"

/******************************************************************
*                        Class Sphere                             *
******************************************************************/
class Sphere: public Intersectable {
public:
    double radius;  // radius of the sphere
    vec3 pos;  // position, emission, color of the sphere
    shared_ptr<Material> mat;
    BoundingBox bbox;

    Sphere(vec3 p, double rad, shared_ptr<Material> mat=nullptr) {
        pos = p;
        radius = rad;
        this->mat = mat;
        bbox = boundingBox();
    }

    void setMaterial(shared_ptr<Material> mat) override {
        this->mat = mat;
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