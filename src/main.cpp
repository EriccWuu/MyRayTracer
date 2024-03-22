#include "iostream"
#include "cmath"
#include <chrono>
#include "Utility.h"
#include "MathLib.h"
#include "Camera.h"

void test() {
    double aspectRatio = 16.0 / 9.0;
    int height = 675, width = (int)(height * aspectRatio);
    TGAImage image(width, height, TGAImage::RGB);
    std::vector<shared_ptr<Intersectable>> world;

    // auto ground_material = make_shared<Lambertian>(vec3(0.5, 0.5, 0.5));
    auto checker = make_shared<SolidColor>(vec3(0.5, 0.5, 0.5));
    auto ground_material = make_shared<Lambertian>(checker);
    world.push_back(make_shared<Sphere>(vec3(0,-1e3,0), 1e3, ground_material));

    for (int a = -11; a < 11; a++) {
        for (int b = -11; b < 11; b++) {
            auto choose_mat = randDouble();
            vec3 center(a + 0.9*randDouble(), 0.2, b + 0.9*randDouble());

            if ((center - vec3(4, 0.2, 0)).norm() > 0.9) {
                shared_ptr<Material> sphere_material;

                if (choose_mat < 0.8) {
                    // diffuse
                    auto albedo = randVec3();
                    sphere_material = make_shared<Lambertian>(albedo);
                    world.push_back(make_shared<Sphere>(center, 0.2, sphere_material));
                } else if (choose_mat < 0.95) {
                    // metal
                    auto albedo = randVec3(0.5, 1);
                    auto fuzz = randDouble(0, 0.5);
                    sphere_material = make_shared<Metal>(albedo, fuzz);
                    world.push_back(make_shared<Sphere>(center, 0.2, sphere_material));
                } else {
                    // glass
                    sphere_material = make_shared<Dielectric>(vec3(1, 1, 1), 1.0 / 1.5);
                    world.push_back(make_shared<Sphere>(center, 0.2, sphere_material));
                }
            }
        }
    }

    auto material1 = std::make_shared<Dielectric>(vec3(1, 1, 1), 1.0 / 1.5);
    world.push_back(std::make_shared<Sphere>(vec3(0, 1, 0), 1.0, material1));

    auto material2 = std::make_shared<Lambertian>(vec3(0.4, 0.2, 0.1));
    world.push_back(std::make_shared<Sphere>(vec3(-4, 1, 0), 1.0, material2));

    auto material3 = std::make_shared<Metal>(vec3(0.7, 0.6, 0.5), 0.0);
    world.push_back(std::make_shared<Sphere>(vec3(4, 1, 0), 1.0, material3));

    double fov     = 30;
    double near = -0.5, far = -50;
    vec3 position  = vec3(13,2,3);
    vec3 direction = ZERO_VEC3 - position;
    vec3 up        = vec3(0,1,0);

    Camera cam(position, direction, up, height, aspectRatio, fov, near, far);

    cam.spp         = 500;
    cam.maxDepth    = 15;
    cam.defocusAngle = 0;
    cam.near    = -0.5;
    cam.backgroundLight = ONE_VEC3 * 0.9;

    std::sort(world.begin(), world.begin() + world.size(), [&](const shared_ptr<Intersectable> &a, const shared_ptr<Intersectable> &b) {return compare(a, b);});

    BVHNode world_(world, 0, world.size());
    cam.render(world_, image);

    image.write_tga_file("result.tga");
}

void test1() {
    double aspectRatio = 16.0 / 9.0;
    int height = 450, width = (int)(height * aspectRatio);
    double fov = 30;
    double near = -1/std::tan(PI/6), far = -50;
    double a = 20;
    vec3 up = {0, 1, 0};
    vec3 camPos = {0, a, 150};
    // vec3 camDir = vec3(0, 0, -1);
    vec3 camDir = vec3(0, a, 0) - camPos;

    TGAImage image(width, height, TGAImage::RGB);
    Camera camera(camPos, camDir, up, height, aspectRatio, fov, near, far);
    camera.spp = 100;
    camera.maxDepth = 10;
    camera.defocusAngle = 0;
    camera.backgroundLight = ONE_VEC3 * 0.8;

    auto matGround = std::make_shared<Lambertian>(vec3(0.8, 0.8, 0.8));
    auto matCeil = std::make_shared<Lambertian>(vec3(0.8, 0.8, 0.8));
    auto matBackWall = std::make_shared<Lambertian>(vec3(0.8, 0.8, 0.8));
    auto matFrontWall = std::make_shared<Lambertian>(vec3(0, 0, 0));
    auto matLeftWall = std::make_shared<Lambertian>(vec3(0.75, 0.25, 0.25));
    auto matRightWall = std::make_shared<Lambertian>(vec3(0.25, 0.25, 0.75));
    auto matCenter = std::make_shared<Lambertian>(vec3(0.7, 0.3, 0.3));
    auto matLight = std::make_shared<Emission>(12*ONE_VEC3);
    auto matGrass = std::make_shared<Dielectric>(vec3(1.0, 1.0, 1.0), 1.0/1.5);
    auto matMirror  = std::make_shared<Metal>(vec3(0.8, 0.8, 0.8), 0);
    auto matDiffuse = std::make_shared<Metal>(vec3(0.8, 0.6, 0.2), 0.3);

    Sphere ground(vec3(0, -1e5, 0), 1e5, matGround);
    Sphere ceil(vec3(0, 1e5 + 2*a, 0), 1e5, matCeil);
    Sphere backWall(vec3(0, 0, -1e5 - a), 1e5, matBackWall);
    Sphere frontWall(vec3(0, 0, 1e5 + 2*a), 1e5, matFrontWall);
    Sphere leftWall(vec3(-1e5 - 2*a, 0, 0), 1e5, matLeftWall);
    Sphere rightWall(vec3(1e5 + 2*a, 0, 0), 1e5, matRightWall);
    Sphere light(vec3(0, 2*a + 140, -a/4), 142, matLight);
    Sphere metalBoll(vec3(-50, 20, -20), 20, matMirror);
    Sphere diffuseBoll(vec3(50, 20, -20), 20, matDiffuse);
    Sphere grassBoll(vec3(0, 20, -20), 20, matGrass);
    Sphere right1(vec3(50, 5, -40), 5, matDiffuse);

    std::vector<shared_ptr<Intersectable>> objects;
    objects.push_back(make_shared<Sphere>(ground));
    // objects.push_back(make_shared<Sphere>(ceil));
    // objects.push_back(make_shared<Sphere>(backWall));
    // objects.push_back(make_shared<Sphere>(frontWall));
    // objects.push_back(make_shared<Sphere>(leftWall));
    // objects.push_back(make_shared<Sphere>(rightWall));
    // objects.push_back(make_shared<Sphere>(light));
    objects.push_back(make_shared<Sphere>(grassBoll));
    objects.push_back(make_shared<Sphere>(metalBoll));
    objects.push_back(make_shared<Sphere>(diffuseBoll));
    // objects.push_back(make_shared<Sphere>(right1));

    std::sort(objects.begin(), objects.begin() + objects.size(), [&](const shared_ptr<Intersectable> &a, const shared_ptr<Intersectable> &b) {return compare(a, b);});

    BVHNode world(objects, 0, objects.size());
    camera.render(world, image);

    image.write_tga_file("result.tga");
}

void earth() {
    double aspectRatio = 16.0 / 9.0;
    int height = 675, width = (int)(height * aspectRatio);
    TGAImage image(width, height, TGAImage::RGB);
    TGAImage earthTexture;
    earthTexture.read_tga_file("assets/earthmap.tga");

    std::vector<shared_ptr<Intersectable>> world;

    auto earthTex = make_shared<ImageTexture>(earthTexture);
    auto earth_material = make_shared<Lambertian>(earthTex);
    auto noiseTex = make_shared<NoiseTexture>(4);
    auto ground_material = make_shared<Lambertian>(noiseTex);
    world.push_back(make_shared<Sphere>(vec3(0, 3, 0), 3, earth_material));
    world.push_back(make_shared<Sphere>(vec3(0, -1e3, 0), 1e3, ground_material));

    double fov     = 60;
    double near = -0.5, far = -50;
    vec3 position  = vec3(0, 3, 10);
    vec3 direction = vec3(0, 2, 0) - position;
    vec3 up        = vec3(0,1,0);

    Camera cam(position, direction, up, height, aspectRatio, fov, near, far);

    cam.spp         = 10;
    cam.maxDepth    = 10;
    cam.defocusAngle = 0;
    cam.near    = -0.5;
    cam.backgroundLight = ONE_VEC3*0.8;

    std::sort(world.begin(), world.begin() + world.size(), [&](const shared_ptr<Intersectable> &a, const shared_ptr<Intersectable> &b) {return compare(a, b);});

    BVHNode world_(world, 0, world.size());
    cam.render(world_, image);

    image.write_tga_file("result.tga");
}

void cornellBox() {
    double aspectRatio = 1.0;
    int height = 600, width = (int)(height * aspectRatio);
    TGAImage image(width, height, TGAImage::RGB);

    vector<shared_ptr<Intersectable>> world;

    TGAImage earthImage;
    earthImage.read_tga_file("assets/earthmap.tga");
    auto earth_material = make_shared<Lambertian>(make_shared<ImageTexture>(earthImage));
    auto red_material = make_shared<Lambertian>(vec3(0.65, 0.05, 0.05));
    auto white_material = make_shared<Lambertian>(vec3(0.73, 0.73, 0.73));
    auto green_material = make_shared<Lambertian>(vec3(0.12, 0.45, 0.15));
    auto light_material = make_shared<Emission>(ONE_VEC3 * 15);
    auto glass_material = make_shared<Dielectric>(ONE_VEC3, 1.0 / 1.5);
    auto sphere_material = make_shared<Lambertian>(vec3(1.0, 0.2, 0.2));
    auto metal_material = make_shared<Metal>(vec3(0.8, 0.6, 0.2), 0.3);
    auto mirror_material = make_shared<Metal>(vec3(0.8, 0.8, 0.8), 0.0);

    double a = 6;
    vec3 Q;
    Q = {-1, a - 1e-3, -1}; // Light
    vector<Vertex> light = {Vertex(Q), Vertex(Q + 2*XA), Vertex(Q + 2*ZA), Vertex(Q + 2*(XA + ZA))};
    vector<int> lightVertInd = {0, 1, 2, 1, 3, 2};
    Q = {-0.5*a, 0, -0.5*a};
    // ground vertics
    vector<Vertex> ground = {Vertex(Q), Vertex(Q + a*ZA), Vertex(Q + a*XA), Vertex(Q + a*(XA + ZA))};
    vector<int> groundVertInd = {0, 1, 2, 1, 3, 2};
    // back-wall vertics
    vector<Vertex> backWall = {Vertex(Q), Vertex(Q + a*XA), Vertex(Q + a*YA), Vertex(Q + a*(XA + YA))};
    vector<int> backWallVertInd = {0, 1, 2, 1, 3, 2};
    // left-wall vertics
    vector<Vertex> leftWall = {Vertex(Q), Vertex(Q + a*YA), Vertex(Q + a*ZA), Vertex(Q + a*(YA + ZA))};
    vector<int> leftWallVertInd = {0, 1, 2, 1, 3, 2};
    Q = {0.5*a, a, 0.5*a};
    // right-wall vertics
    vector<Vertex> rightWall = {Vertex(Q), Vertex(Q - a*ZA), Vertex(Q - a*YA), Vertex(Q - a*(YA + ZA))};
    vector<int> rightWallVertInd = {0, 1, 2, 1, 3, 2};
    // floor vertics
    vector<Vertex> floor = {Vertex(Q), Vertex(Q - a*XA), Vertex(Q - a*ZA), Vertex(Q - a*(XA + ZA))};
    vector<int> floorVertInd = {0, 1, 2, 1, 3, 2};
    // Mirror vertics
    // Q = vec3(0.5*a - 1e-3, a, 0.5*a) + vec3(0.5*a - 1e-3, -1, -1);
    Q = vec3(0.5*a - 1e-3, a, 0.5*a) +  vec3(0, -1, -1);
    vector<Vertex> mirror = {Vertex(Q), Vertex(Q - 4*ZA), Vertex(Q - 4*YA), Vertex(Q - 4*(YA + ZA))};
    vector<int> mirrorVertInd = {0, 1, 2, 1, 3, 2};
    Mesh box1 = rectangle(1.8, 3.5);
    box1.setMaterial(white_material);
    Mesh box2 = rectangle(1.5, 1.5);
    box2.setMaterial(white_material);

    Sphere metalBall(ZERO_VEC3, 0.5, metal_material);
    Sphere earthBall(ZERO_VEC3, 0.8, earth_material);
    Sphere glassBall(ZERO_VEC3, 0.8, glass_material);

    mat4 S = scale(1.0);
    mat4 R = rotateY(20);
    mat4 T = translate(vec3(-1, 0, -1));
    mat4 SRT = T*R*S;
    box1.transform(SRT);

    R = rotateY(-30);
    T = translate(vec3(1, 0, 1));
    SRT = T*R*S;
    box2.transform(SRT);

    T = translate(vec3(1, 2, 1));
    metalBall.pos = (T*vec4(metalBall.pos, 1.0)).xyz();

    T = translate(vec3(1.8, 4, -1.5));
    earthBall.pos = (T*vec4(earthBall.pos, 1.0)).xyz();

    T = translate(vec3(-1, 1.5, 1.5));
    glassBall.pos = (T*vec4(glassBall.pos, 1.0)).xyz();

    world.push_back(make_shared<Mesh>(light, lightVertInd, light_material));
    world.push_back(make_shared<Mesh>(ground, groundVertInd, white_material));
    world.push_back(make_shared<Mesh>(backWall, backWallVertInd, white_material));
    world.push_back(make_shared<Mesh>(leftWall, leftWallVertInd, green_material));
    world.push_back(make_shared<Mesh>(rightWall, rightWallVertInd, red_material));
    world.push_back(make_shared<Mesh>(floor, floorVertInd, white_material));
    world.push_back(make_shared<Mesh>(mirror, mirrorVertInd, mirror_material));
    world.push_back(make_shared<Mesh>(box1));
    world.push_back(make_shared<Mesh>(box2));
    world.push_back(make_shared<Sphere>(metalBall));
    world.push_back(make_shared<Sphere>(earthBall));
    world.push_back(make_shared<Sphere>(glassBall));

    // vec3 in = vec3(0.9, -0.1, 0).normalize();
    // vec3 n = vec3(0, 1, 0);
    // vec3 r = refract(in, n, 1.0 / 1.5);
    // double sinTheta1 = sqrt(1 - (in*n)*(in*n));
    // double sinTheta2 = sqrt(1 - (r*n)*(r*n));
    // std::cout << sinTheta1 << ' ' << sinTheta2*1.5 << '\n';

    double fov     = 60;
    double near = -0.5, far = -50;
    vec3 position  = vec3(0, 0.5*a, 8);
    vec3 direction = vec3(0, 0.5*a, 0) - position;
    vec3 up        = vec3(0, 1, 0);

    Camera cam(position, direction, up, height, aspectRatio, fov, near, far);

    cam.spp         = 500;
    cam.maxDepth    = 10;
    cam.defocusAngle = 0;
    cam.backgroundLight = ZERO_VEC3;
    cam.SSAA = 2;

    std::sort(world.begin(), world.begin() + world.size(), [&](const shared_ptr<Intersectable> &a, const shared_ptr<Intersectable> &b) {return compare(a, b);});

    BVHNode world_(world, 0, world.size());
    cam.render(world_, image);

    image.write_tga_file("result.tga");
}

int main() {
    int idx = 3;
    auto start_time = std::chrono::high_resolution_clock::now();
    switch (idx) {
        case 0: test(); break;
        case 1: test1(); break;
        case 2: earth(); break;
        case 3: cornellBox(); break;
    }
    auto end_time = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(end_time - start_time);
    std::cout << "\nTime cost: " << duration.count() << " ms" << std::endl;

    return 0;
}