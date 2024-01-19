#ifndef __CAMERA_H__
#define __CAMERA_H__

#include "MathLib.h"
#include "TGAImage.h"
#include "Utility.h"

class Camera {
public:
    vec4 position;
    vec4 direction;
    vec4 up;
    vec4 right;
    double fov;
    double aspectRatio;
    double near, far;
    int screen_w, screen_h;

    TGAImage frame;
    TGAImage zbuffer;

    Camera() = default;
    Camera(vec3 pos, vec3 dir, vec3 u, int h=100, double r=1, double fov=60, double n=-0.5, double f=-50);

    void set(vec3 pos, vec3 dir, vec3 u, double h=100, double r=1, double fov=60, double n=-0.5, double f=-50);
    mat4 view();
    mat4 projection();
    mat4 viewport();
    void init();
    void render(Interlist &objects, TGAImage &image);

private:
    vec3 view_o;
    vec3 viewport_o, viewport_u, viewport_v;

    inline double width();
    inline double height();
    Ray  getRay(const int &i, const int &j);
    vec3 pixSampleSquare();
    vec3 pixSampleDisk(double radius);
    bool intersect_all(Interlist &objects, const Ray &ray, Interval rayt, Inter_record &rec);
};

#endif