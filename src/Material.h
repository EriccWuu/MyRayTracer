#ifndef __MATERIAL_H__
#define __MATERIAL_H__

#include "MathLib.h"
#include "Texture.h"
#include "Ray.h"

class Ray;
struct InterRecord;

// material types, used in radiance()
enum ReflType { 
    DIFF, SPEC, REFR
};

class Material {
public:
    virtual ~Material() = default;
    virtual bool scatter(
        const Ray &rayIn,
        const InterRecord &rec,
        Ray &scattered,
        vec3 &attenuation
    ) const = 0;
    virtual vec3 emit() const {
        return ZERO_VEC3;
    }
    virtual vec3 color() const = 0;
};

class Lambertian: public Material {
public:
    Lambertian(const vec3 &a): albedo(std::make_shared<SolidColor>(a)) {}
    Lambertian(shared_ptr<Texture> a): albedo(a) {}
    bool scatter(const Ray &rayIn, const InterRecord &rec, Ray &scattered, vec3 &attenuation) const override {
        vec3 n = rec.inward ? rec.normal : -rec.normal;
        vec3 scatterDir = randVecSemisphere(n);    // Lambertian Reflection
        scattered = Ray(rec.p, scatterDir);
        attenuation = albedo->value(rec.uv.x, rec.uv.y, rec.p);
        return true;
    }
    vec3 color() const override { return albedo->value(0, 0, vec3(1, 0, 0)); }

private:
    shared_ptr<Texture> albedo;
};

class Metal: public Material {
public:
    Metal(const vec3 &a, double f): albedo(a), fuzzy(clamp(f)) {}
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
        double refractRatio = rec.inward ? ir : 1/ir;
        vec3 n = rec.inward ? rec.normal : -rec.normal;
        vec3 dir;
        dir = refract(rayin.dir, n, refractRatio); // Refract
        scattered = Ray(rec.p, dir);
        attenuation = albedo;
        return true;
    }
    vec3 color() const override { return albedo; }

private:
    vec3 albedo;
    double ir;
};

class Emission: public Material {
public:
    Emission(const vec3 &e):emission(e) {}
    bool scatter(const Ray &rayIn, const InterRecord &rec, Ray &scattered, vec3 &attenuation) const override {
        return false;
    }
    vec3 emit() const override {
        return emission;
    }
    vec3 color() const override { return emission; }

private:
    vec3 emission;
};
#endif