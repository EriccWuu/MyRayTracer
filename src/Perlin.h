#ifndef __PERLIN_H__
#define __PERLIN_H__

#include "MathLib.h"

class Perlin {
public:
    Perlin() {
        randvec = new vec3[cnt];
        for (int i = 0; i < cnt; i ++) {
            randvec[i] = randVec3(-1, 1);
        }
        permX = perlinGeneratePerm();
        permY = perlinGeneratePerm();
        permZ = perlinGeneratePerm();
    }

    ~Perlin() {
        delete[] randvec;
        delete[] permX;
        delete[] permY;
        delete[] permZ;
    }

    double noise(const vec3 & p) const {
        auto u = p.x - floor(p.x);
        auto v = p.y - floor(p.y);
        auto w = p.z - floor(p.z);
        // Hermitian Smoothing
        u = u*u*(3-2*u);
        v = v*v*(3-2*v);
        w = w*w*(3-2*w);
        auto i = static_cast<int>(floor(p.x));
        auto j = static_cast<int>(floor(p.y));
        auto k = static_cast<int>(floor(p.z));
        vec3 c[2][2][2];

        for (int ii = 0; ii < 2; ii ++)
            for (int jj = 0; jj < 2; jj ++)
                for (int kk = 0; kk < 2; kk ++)
                    c[ii][jj][kk] = randvec[permX[(ii+i)&255] ^ permY[(jj+j)&255] ^ permZ[(kk+k)&255]];

        return trilinearInterp(c, u , v, w);
    }

    double turbulance(const vec3 &p, int depth = 7) const {
        auto accum = 0.;
        auto tmp = p;
        auto weight = 1.;

        for (int i = 0; i < depth; i ++) {
            accum += weight * noise(tmp);
            weight *= 0.5;
            tmp *=2;
        }

        return fabs(accum);
    }

private:
    static const int cnt = 256;
    vec3* randvec;
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

    static double trilinearInterp(vec3 c[2][2][2], double u, double v, double w) {
        double accum = 0.0;
        for (int i = 0; i < 2; i ++)
            for (int j = 0; j < 2; j ++)
                for (int k = 0; k < 2; k ++) {
                    vec3 weight(u-i, v-j, w-k);
                    accum += (i*u + (1-i)*(1-u))*(j*v + (1-j)*(1-v))*(k*w + (1-k)*(1-w))*c[i][j][k]*weight;
                }
        return accum;
    }

};

#endif