#include "iostream"
#include "Matrix.h"
#include "Camera.h"
#include "utility.h"

vec3 origin = {0, 0, 0};
double aspectRatio = 16.0 / 9.0;
int width = 800, height = (int)(width / aspectRatio);
double view_aspectRatio = (double)width / height;
double fov = 60;
double near = -1, far = -50;
vec3 up = {0, 0, 1};
vec3 camPos = {0, 0, 3};
vec3 camDir = origin - camPos;

int main() {
    TGAImage image(width, height, TGAImage::RGB);
    Camera Camera(camPos, camDir, up, width, height, fov, view_aspectRatio, near, far);
    vec3 v = {1, 2, 3};
    vec3 u = {1, 2, 1};
    mat3 M = {{{1, 2, 3} ,{4, 5, 6}, {7, 8, 8}}};
    std::cout << cross(u, v)  << '\n';
    std::cout << M.invert() << '\n';

    std::cout << width << ' ' << height << '\n';
    std::cout << aspectRatio << ' ' << view_aspectRatio << '\n';
    
    return 0;
}