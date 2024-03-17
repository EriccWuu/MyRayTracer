#ifndef __TEXTURE_H__
#define __TEXTURE_H__

#include "MathLib.h"
#include "TGAImage.h"
#include "Perlin.h"

class Texture {
public:
    virtual ~Texture() = default;
    virtual vec3 value(double u, double v, const vec3 &p) const = 0;
};

class SolidColor: public Texture {
public:
    SolidColor(vec3 c): color(c) {}
    SolidColor(double r, double g, double b): color(vec3(r, g, b)) {}
    vec3 value(double u, double v, const vec3 &p) const override {
        return color;
    }

private:
    vec3 color;
};

class ImageTexture: public Texture {
public:
    ImageTexture(const std::string filepath) {
        texture.read_tga_file(filepath);
    }
    ImageTexture(TGAImage tex): texture(tex) {}

    vec3 value(double u, double v, const vec3 &p) const override {
        int i = static_cast<int>(u * texture.width());
        int j = texture.height() - static_cast<int>(v * texture.height());
        TGAColor c = texture.get(i, j);
        return vec3(c[2], c[1], c[0]) / 255.0;
    }

private:
    TGAImage texture;
};

class NoiseTexture: public Texture {
public:
    NoiseTexture() {}
    NoiseTexture(double sc): scale(sc) {}

    vec3 value(double u, double v, const vec3 &p) const override {
        auto s = scale * p;
        // return noise.turbulance(s) * ONE_VEC3;
        return 0.5 * (1 + sin(s.z + 10*noise.turbulance(s))) * ONE_VEC3;
    }

private:
    Perlin noise;
    double scale;   // to tweaking the frequency
};

#endif