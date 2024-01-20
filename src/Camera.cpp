#include "Camera.h"
#include "stdio.h"

/******************************************************************
*                    Public Method Defination                     *
******************************************************************/
Camera::Camera(vec3 pos, vec3 dir, vec3 u, int h, double r, double fov, double n, double f) {
    position = vec4(pos, 1.f);
    direction = vec4(dir.normalize());
    right = vec4(cross(direction.xyz(), u).normalize());
    up = vec4(cross(right.xyz(), direction.xyz()).normalize());
    this->fov = fov;
    screen_h = h;
    screen_w = screen_h * r;
    aspectRatio = double(screen_w) / screen_h;
    near = n;
    far = f;
}

void Camera::set(vec3 pos, vec3 dir, vec3 u, double h, double r, double fov, double n, double f) {
    position = vec4(pos, 1.f);
    direction = vec4(dir.normalize());
    right = vec4(cross(direction.xyz(), u).normalize());
    up = vec4(cross(right.xyz(), direction.xyz()).normalize());
    this->fov = fov;
    screen_h = h;
    screen_w = screen_h * r;
    aspectRatio = double(screen_w) / screen_h;
    near = n;
    far = f;
}

mat4 Camera::view() {
    vec3 z = -direction.xyz().normalize();
    vec3 x = cross(up.xyz(), z).normalize();
    vec3 y = cross(z, x).normalize();

    mat4 Rview = E44;
    mat4 Tview = E44;
    Tview[0][3] = -position.x;
    Tview[1][3] = -position.y;
    Tview[2][3] = -position.z;
    Rview[0] = vec4(x);
    Rview[1] = vec4(y);
    Rview[2] = vec4(z);
    return Rview * Tview;
}

mat4 Camera::projection() {
    double h = -2*near*tan(deg2rad(fov/2));
    double w = aspectRatio * h;
    double n = near;
    double f = far;
    mat4 ortho, persp2ortho;
    ortho = {{{2.0/w, 0.f, 0.f, 0.f},
                    {0.f, 2.0/h, 0.f, 0.f},
                    {0.f, 0.f, 2.0/(n-f), -(n+f)/(n-f)},
                    {0.f, 0.f, 0.f, 1.f}}};
    persp2ortho = {{{n, 0.f, 0.f, 0.f},
                        {0.f, n, 0.f, 0.f},
                        {0.f, 0.f, n+f, -n*f},
                        {0.f, 0.f, 1.f, 0.f}}};
    return ortho * persp2ortho;
}

mat4 Camera::viewport() {
    mat4 viewport = E44;
    viewport[0][3] = screen_w / 2.f;
    viewport[1][3] = screen_h / 2.f;
    viewport[0][0] = screen_w / 2.f;
    viewport[1][1] = screen_h / 2.f;
    return viewport;
}

void Camera::init() {
    viewO = vec3(0, 0, 0);
    viewportO = ZA * near;
    viewportU = XA * width() / screen_w;
    viewportV = -YA * height() / screen_h;
}

void Camera::render(Interlist &objects, TGAImage &image) {
    init();

    for (int j = 0; j < screen_h; j ++) {
        fprintf(stderr, "\rRenderering (%d spp): %5.2f%%", spp, 100.*j/(screen_h-1));
        for (int i = 0; i < screen_w; i ++) {
            vec3 color(0, 0, 0);
            for (int k = 0; k < spp; k ++) {
                Ray ray = getRay(i, j);
                color += radiance(ray, 0, objects);
            }
            color = clampVec3(color / spp); // Clamp to [0, 1]
            gammaCorrection(color); // Gamma correction
            image.set(i, j, TGAColor(255*color));
        }
    }
    std::cout << '\n';
}

/******************************************************************
*                    Private Method Defination                    *
******************************************************************/
inline double Camera::width() {
    return aspectRatio * height();
}

inline double Camera::height() {
    return -2*near*std::tan(deg2rad(fov/2));
}

Ray Camera::getRay(const int &i, const int &j) {
    vec3 pixPos = viewportO + (i - (screen_w >> 1) + 0.5)*viewportU + (j - (screen_h >> 1) + 0.5)*viewportV;
    pixPos += pixSampleSquare();
    return Ray(viewO, pixPos);
}

vec3 Camera::radiance(const Ray &r, int depth, const Interlist &obj) {
    InterRecord rec;

    // Return if no intersection
    if (!intersect_all(obj, r, Interval(0.001, INF), rec)) {
        double a = 0.5*(r.direction().normalize().y + 1.0);
        return ((1.0 - a)*vec3(1.0, 1.0, 1.0) + a*vec3(0.5,0.7,1.0));
    }
    // vec3 c = rec.mat->color();
    // double p = c.x>c.y && c.x>c.z ? c.x : c.y>c.z ? c.y : c.z; // max refl
    if (++depth > maxDepth) {
        return ZERO_VEC3;
    }

    Ray scattered;
    vec3 attenuation;
    rec.mat->scatter(r, rec, scattered, attenuation);
    return radiance(scattered, depth, obj).mult(attenuation);
}

inline vec3 Camera::pixSampleSquare() {
    vec2 sample = sampleUnitSpuare();
    double sx = sample.x - 0.5;
    double sy = sample.y - 0.5;
    return viewportU*sx + viewportV*sy;
}

inline vec3 Camera::pixSampleDisk(double radius) {
    vec2 sample = sampleUnitDisk();
    sample *= radius;
    return viewportU*sample.x + viewportV*sample.y;
}

bool Camera::intersect_all(const Interlist &obj, const Ray &ray, Interval rayt, InterRecord &rec) {
    InterRecord tmp_rec;
    bool hit_anything = false;
    
    for (const auto& obj: obj) {
        if (obj->intersect(ray, rayt, tmp_rec)) {
            hit_anything = true;
            rayt.max = tmp_rec.t;
            rec = tmp_rec;
        }
    }

    return hit_anything;
}

