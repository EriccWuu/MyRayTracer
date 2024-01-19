#include "Camera.h"

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

    mat4 Rview = mat4::identity();
    mat4 Tview = mat4::identity();
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
    mat4 viewport = mat4::identity();
    viewport[0][3] = screen_w / 2.f;
    viewport[1][3] = screen_h / 2.f;
    viewport[0][0] = screen_w / 2.f;
    viewport[1][1] = screen_h / 2.f;
    return viewport;
}

void Camera::init() {
    view_o = vec3(0, 0, 0);
    viewport_o = vec3(0, 0, 1) * near;
    viewport_u = vec3(1, 0, 0) * width() / screen_w;
    viewport_v = vec3(0, 1, 0) * height() / screen_h;
}

void Camera::render(Interlist &objects, TGAImage &image) {
    init();
    Inter_record rec;

    for (int j = 0; j < screen_h; j ++) {
        std::cout << "\rScanlines remaining: " << (screen_h - 1 - j) << ' ' << std::flush;
        for (int i = 0; i < screen_w; i ++) {
            Ray ray = getRay(i, j);
            double a = 0.5*(ray.direction().y + 1.0);
            vec3 color;
            if (!intersect_all(objects, ray, Interval(0, INF), rec)) {
                color = 255*((1.0 - a)*vec3(1.0, 1.0, 1.0) + a*vec3(0.5,0.7,1.0));
                image.set(i, screen_h-j, TGAColor(color));
                continue;
            }
            color = 0.5*(rec.normal + 1);
            image.set(i, screen_h-j, TGAColor(255*color));
        }
    }
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
    vec3 pixPos = viewport_o + (i - (screen_w >> 1) + 0.5)*viewport_u + (j - (screen_h >> 1) + 0.5)*viewport_v;
    pixPos += pixSampleSquare();
    return Ray(view_o, pixPos);
}

inline vec3 Camera::pixSampleSquare() {
    double dx = randDouble() - 0.5;
    double dy = randDouble() - 0.5;
    return viewport_u*dx + viewport_v*dy;
}

// vec3 Camera::pixSampleDisk(double radius) {

// }

bool Camera::intersect_all(Interlist &objects, const Ray &ray, Interval rayt, Inter_record &rec) {
    Inter_record tmp_rec;
    bool hit_anything = false;
    
    for (const auto& obj: objects) {
        if (obj->intersect(ray, rayt, tmp_rec)) {
            hit_anything = true;
            rayt.max = tmp_rec.t;
            rec = tmp_rec;
        }
    }
    return hit_anything;
}

