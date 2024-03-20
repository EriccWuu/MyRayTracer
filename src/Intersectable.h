#ifndef __INTERSECTABLE_H__
#define __INTERSECTABLE_H__

#include "MathLib.h"
#include "Ray.h"
#include "Material.h"

/******************************************************************
*                        Class Boundingbox                        *
******************************************************************/
class BoundingBox {
public:
    vec3 pmin, pmax;

    BoundingBox(): pmin(ZERO_VEC3), pmax(ZERO_VEC3) {}
    BoundingBox(const vec3 &pmin_, const vec3 &pmax_) {
        pmin = pmin_;
        pmax = pmax_;
        pad(pmin, pmax);
    }
    BoundingBox(const BoundingBox &a, const BoundingBox &b) {
        float ax = std::min(a.pmin.x, b.pmin.x);
        float ay = std::min(a.pmin.y, b.pmin.y);
        float az = std::min(a.pmin.z, b.pmin.z);
        float bx = std::max(a.pmax.x, b.pmax.x);
        float by = std::max(a.pmax.y, b.pmax.y);
        float bz = std::max(a.pmax.z, b.pmax.z);
        pmin = vec3(ax, ay, az);
        pmax = vec3(bx, by, bz);
        pad(pmin, pmax);
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

private:
    void pad(vec3 &pmin, vec3 &pmax) {
        auto delta = 1e-3;
        pmax.x = (pmax.x - pmin.x > delta) ? pmax.x : pmax.x + delta;
        pmax.y = (pmax.y - pmin.y > delta) ? pmax.y : pmax.y + delta;
        pmax.z = (pmax.z - pmin.z > delta) ? pmax.z : pmax.z + delta;
    }
};

/******************************************************************
*                        Class Intersectable                      *
******************************************************************/
class Intersectable {
public:
    // BoundingBox bbox;
    virtual ~Intersectable() = default;
    // Compute wheather a ray hitted the object
    virtual bool intersect(const Ray &ray, Interval rayt, InterRecord &rec) const = 0;
    // Compute axis-aligned bounding box of the object
    virtual BoundingBox boundingBox() const = 0;
    // Set material for Intersectable instance
    virtual void setMaterial(std::shared_ptr<Material> mat) {}
    // Update data in Intersectable instance
    virtual void update() {}
};

static bool compare(const shared_ptr<Intersectable> &a, const shared_ptr<Intersectable> &b, int axis = 0) {
    return a->boundingBox().pmin[axis] < b->boundingBox().pmin[axis];
}

#endif