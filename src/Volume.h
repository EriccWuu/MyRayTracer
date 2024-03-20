#ifndef __VOLUME_H__
#define __VOLUME_H__

#include "Intersectable.h"
#include "Triangle.h"

/******************************************************************
*                        Class Volume                             *
******************************************************************/
class Volume: public Intersectable {
public:
    std::vector<Vertex> *verts;
    int vertIndex[3];
    shared_ptr<Material> mat;

    Volume(shared_ptr<Intersectable> b, double d, shared_ptr<Texture> tex)
    : boundary(b), neg_inv_density(-1.0/d), phaseFunction(std::make_shared<Isotropic>(tex)) {}

    Volume(shared_ptr<Intersectable> b, double d, const vec3 tex)
    : boundary(b), neg_inv_density(-1.0/d), phaseFunction(std::make_shared<Isotropic>(tex)) {}


    void setMaterial(shared_ptr<Material> mat) {
        this->mat = mat;
    }

    BoundingBox boundingBox() const override {
        return boundary->boundingBox();
    }

    vec2 uv(const vec3 &p) {
        return vec2(0, 0);
    }

    // Compute wheather a ray hitted the Triangle
    bool intersect(const Ray &ray, Interval rayt, InterRecord &rec) const override {
        InterRecord rec1, rec2;
        if (!boundary->intersect(ray, Interval::universe, rec1))
            return false;
        if (!boundary->intersect(ray, Interval(rec1.t + 1e-3, INF), rec2))
            return false;

        if (rec1.t < rayt.min) rec1.t = rayt.min;
        if (rec2.t > rayt.max) rec2.t = rayt.max;

        if (rec1.t >= rec2.t)
            return false;

        if (rec1.t < 0)
            rec1.t = 0;

        double rayLength = ray.dir.norm();
        double distInsideBoundry = (rec2.t - rec1.t) * rayLength;
        double hitDistance = neg_inv_density * log(randDouble());

        if (hitDistance > distInsideBoundry)
            return false;
        rec.t = rec1.t + hitDistance / rayLength;
        rec.p = ray.at(rec.t);
        rec.normal = ONE_VEC3;
        rec.inward = true;
        rec.mat = phaseFunction;

        return true;
    }

private:
    shared_ptr<Intersectable> boundary;
    double neg_inv_density;
    shared_ptr<Isotropic> phaseFunction;

};
#endif