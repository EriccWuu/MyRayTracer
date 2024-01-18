#include "iostream"
#include "cmath"
#include "Matrix.h"
#include "Camera.h"
#include "utility.h"

vec3 origin = {0, 0, 0};
double aspectRatio = 16.0 / 9.0;
int width = 800, height = (int)(width / aspectRatio);
double view_aspectRatio = (double)width / height;
double fov = 60;
double near = -1/std::tan(PI/6), far = -50;
vec3 up = {0, 1, 0};
vec3 camPos = {0, 0, 3};
vec3 camDir = origin - camPos;

Sphere sphere(2, vec3(0, 0, -5), vec3(1, 0, 0));

int main() {
    TGAImage image(width, height, TGAImage::RGB);
    Camera Camera(camPos, camDir, up, width, height, fov, view_aspectRatio, near, far);

    std::cout << Camera.near << '\n';
    std::cout << Camera.width() << ' ' << Camera.height() << '\n';
    std::cout << Camera.viewport_u() << ' ' << Camera.viewport_v() << '\n';

    vec3 vp_o = vec3(0, 0, 1) * near;
    vec3 vp_u = Camera.viewport_u(), vp_v = Camera.viewport_v();
    for (int j = 0; j < height; j ++) {
        std::cout << "\rScanlines remaining: " << (height - j) << ' ' << std::flush;
        for (int i = 0; i < width; i ++) {
            vec3 rayDir = (vp_o + (i - width/2)*vp_u + (j - height/2)*vp_v).normalize();
            Ray ray(vec3(0,0,0), rayDir);
            double a = 0.5*(ray.direction().y + 1.0);
            vec3 color = 255*((1.0 - a)*vec3(1.0, 1.0, 1.0) + a*vec3(0.5,0.7,1.0));
            if (!sphere.intersect(ray)) {
                image.set(i, height-j, TGAColor(color));
                continue;
            }
            image.set(i, height-j, TGAColor(255*sphere.color));
        }
    }

    image.write_tga_file("result.tga");

    return 0;
}