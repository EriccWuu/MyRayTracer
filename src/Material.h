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
    virtual vec3 color() const = 0;
};

class Lambertian: public Material {
public:
    // Lambertian(const vec3 &a): albedo(a) {}
    Lambertian(const vec3 &a) { albedo = a; }
    bool scatter(const Ray &rayIn, const InterRecord &rec, Ray &scattered, vec3 &attenuation) const override {
        vec3 n = rec.inward ? rec.normal : -rec.normal;
        // vec3 scatterDir = randVecSemisphere(n);    // Random Reflection
        vec3 scatterDir = n + randVecSphere();  // Lambertian Reflection
        auto norm = scatterDir.norm();
        if (norm < EPS) scatterDir = rec.normal;
        else scatterDir /= norm;
        scattered = Ray(rec.p, scatterDir);
        attenuation = albedo;
        return true;
    }
    vec3 color() const override { return albedo; }

private:
    vec3 albedo;
};

class Matel: public Material {
public:
    // matel(const vec3 &a, double f): albedo(a), fuzzy(clamp(f)) {}
    Matel(const vec3 &a, double f) { albedo = a, fuzzy = clamp(f); }
    bool scatter(const Ray &rayin, const InterRecord &rec, Ray &scattered, vec3 &attenuation) const override {
        vec3 n = rec.inward ? rec.normal : -rec.normal;
        vec3 reflectdir = reflect(rayin.dir, n) + fuzzy*randVecSemisphere(n); // reflect direction
        scattered = Ray(rec.p, reflectdir.normalize());
        attenuation = albedo;
        return true;
    }
    vec3 color() const override { return albedo; }

private:
    vec3 albedo;
    double fuzzy;
};

class Dielectric: public Material {
public:
    Dielectric(const vec3 &a, double ir): albedo(a), ir(ir)  {}
    bool scatter(const Ray &rayin, const InterRecord &rec, Ray &scattered, vec3 &attenuation) const override {

        double refractRatio = rec.inward ? 1/ir : ir;
        vec3 n = rec.inward ? rec.normal : -rec.normal;
        double cos_theta = fmin((-rayin.dir * n) , 1.0);
        double sin_theta = sqrt(1.0 - cos_theta*cos_theta);

        vec3 dir;
        dir = refract(rayin.dir, n, refractRatio); // Refract
        // if (refractRatio * sin_theta > 1.0)
        //     dir = reflect(rayin.dir, n);    // Total Reflect
        // else
        //     dir = refract(rayin.dir, n, refractRatio); // Refract
        scattered = Ray(rec.p, dir);
        attenuation = albedo;
        return true;
    }
    vec3 color() const override { return albedo; }

private:
    vec3 albedo;
    double ir;
};
#endif