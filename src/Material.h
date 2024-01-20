#ifndef __MATERIAL_H__
#define __MATERIAL_H__

#include "MathLib.h"
#include "Ray.h"

class Ray;
struct InterRecord;

// material types, used in radiance()
enum ReflType { 
    DIFF, SPEC, REFR
};

class Material {
public:
    // vec3 albedo;

    virtual ~Material() = default;
    virtual bool scatter(
        const Ray &rayIn,
        const InterRecord &rec,
        Ray &scattered,
        vec3 &attenuation
    ) const = 0;
    // vec3 color() { return albedo; }
};

class Lambertian: public Material {
public:
    // Lambertian(const vec3 &a): albedo(a) {}
    Lambertian(const vec3 &a) { albedo = a; }
    bool scatter(const Ray &rayIn, const InterRecord &rec, Ray &scattered, vec3 &attenuation) const override {
        vec3 n = rec.inward ? rec.normal : -rec.normal;
        vec3 scatterDir = randVecSemisphere(n);    // Random Reflection
        // vec3 scatterDir = n + randVecSphere();  // Lambertian Reflection
        if (scatterDir.norm() < EPS) scatterDir = rec.normal;
        scattered = Ray(rec.p, scatterDir);
        attenuation = albedo;
        return true;
    }
    vec3 color() { return albedo; }

private:
    vec3 albedo;
};

class Matel: public Material {
public:
    // Matel(const vec3 &a, double f): albedo(a), fuzzy(clamp(f)) {}
    Matel(const vec3 &a, double f) { albedo = a, fuzzy = clamp(f); }
    bool scatter(const Ray &rayIn, const InterRecord &rec, Ray &scattered, vec3 &attenuation) const override {
        vec3 n = rec.inward ? rec.normal : -rec.normal;
        vec3 reflectDir = reflect(rayIn.dir, n) + fuzzy*randVecSemisphere(n);    // Mirror Reflection
        scattered = Ray(rec.p, reflectDir);
        attenuation = albedo;
        return true;
    }
    vec3 color() { return albedo; }

private:
    vec3 albedo;
    double fuzzy;
};

#endif