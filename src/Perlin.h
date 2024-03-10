#ifndef __PERLIN_H__
#define __PERLIN_H__

#include "MathLib.h"

class Perlin {
public:
    Perlin() {
        randf = new double[cnt];
        for (int i = 0; i < cnt; i ++) {
            randf[i] = randDouble();
        }
        permX = perlinGeneratePerm();
        permY = perlinGeneratePerm();
        permZ = perlinGeneratePerm();
    }

    ~Perlin() {
        delete[] randf;
        delete[] permX;
        delete[] permY;
        delete[] permZ;
    }

    double noise(const vec3 & p) const {
        int i = static_cast<int>(4*p.x) & 255;
        int j = static_cast<int>(4*p.y) & 255;
        int k = static_cast<int>(4*p.z) & 255;

        return randf[permX[i] ^ permY[j] ^ permZ[k]];
    }

private:
    static const int cnt = 256;
    double* randf;
    int* permX;
    int* permY;
    int* permZ;

    static int* perlinGeneratePerm() {
        auto p = new int[cnt];
        for (int i = 0; i < Perlin::cnt; i ++) p[i] = i;

        permute(p, cnt);

        return p;
    }

    static void permute(int *p, int n) {
        for (int i = 0; i < n; i ++) {
            int target = randInt(0, i);
            std::swap(p[target], p[i]);
        }
    }

};

#endif