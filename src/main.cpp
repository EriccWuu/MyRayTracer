#include "iostream"
#include "cmath"
#include <chrono>
#include "MathLib.h"
#include "Camera.h"

vec3 origin = {0, 0, 0};
double aspectRatio = 16.0 / 9.0;
int height = 450, width = (int)(height * aspectRatio);
double fov = 60;
double near = -1/std::tan(PI/6), far = -50;
vec3 up = {0, 1, 0};
vec3 camPos = {0, 0, 3};
vec3 camDir = origin - camPos;

auto matGround = std::make_shared<Lambertian>(vec3(0.8, 0.8, 0.8));
auto matLeftWall = std::make_shared<Lambertian>(vec3(0.75, 0.25, 0.25));
auto matRightWall = std::make_shared<Lambertian>(vec3(0.25, 0.25, 0.75));
auto matCenter = std::make_shared<Lambertian>(vec3(0.7, 0.3, 0.3));
auto matLeft  = std::make_shared<Matel>(vec3(0.8, 0.8, 0.8), 0);
auto matRight = std::make_shared<Matel>(vec3(0.8, 0.6, 0.2), 0.3);

Sphere ground(vec3(0, -1e5 -2, -5), 1e5, matGround);
Sphere leftWall(vec3(-1e5 - 15, 0, -5), 1e5, matLeftWall);
Sphere rightWall(vec3(1e5 + 15, 0, -5), 1e5, matRightWall);
Sphere sphere(vec3(0, 0, -5), 2, matCenter);
Sphere left(vec3(-4, 0, -5), 1.5, matLeft);
Sphere right(vec3(4, 0, -5), 1.5, matRight);
Interlist objects;

int main() {
    TGAImage image(width, height, TGAImage::RGB);
    Camera camera(camPos, camDir, up, height, aspectRatio, fov, near, far);
    camera.spp = 100;
    camera.maxDepth = 5;
    objects.push_back(make_shared<Sphere>(ground));
    // objects.push_back(make_shared<Sphere>(leftWall));
    // objects.push_back(make_shared<Sphere>(rightWall));
    objects.push_back(make_shared<Sphere>(sphere));
    objects.push_back(make_shared<Sphere>(left));
    objects.push_back(make_shared<Sphere>(right));

    auto start_time = std::chrono::high_resolution_clock::now();

    camera.render(objects, image);

    auto end_time = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(end_time - start_time);
    std::cout << "Time cost: " << duration.count() << " ms" << std::endl;

    image.write_tga_file("result.tga");

    return 0;
}