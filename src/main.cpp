#include "iostream"
#include "cmath"
#include "MathLib.h"
#include "Camera.h"
#include "Utility.h"

vec3 origin = {0, 0, 0};
double aspectRatio = 16.0 / 9.0;
int height = 450, width = (int)(height * aspectRatio);
double fov = 60;
double near = -1/std::tan(PI/6), far = -50;
vec3 up = {0, 1, 0};
vec3 camPos = {0, 0, 3};
vec3 camDir = origin - camPos;

Sphere sphere(vec3(0, 0, -5), 2);
Sphere ground(vec3(0, -100, -25), 100);
Interlist objects;

int main() {
    TGAImage image(width, height, TGAImage::RGB);
    Camera camera(camPos, camDir, up, height, aspectRatio, fov, near, far);
    objects.push_back(make_shared<Sphere>(sphere));
    objects.push_back(make_shared<Sphere>(ground));

    camera.render(objects, image);

    image.write_tga_file("result.tga");

    return 0;
}