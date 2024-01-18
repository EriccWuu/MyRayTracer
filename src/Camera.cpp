#include "Camera.h"

Camera::Camera() {
    position = vec4(0.f, 0.f, 0.f, 1.f);
    direction = vec4(0.f, 0.f, -1.f, 0.f);
    up = vec4(0.f, 1.f, 0.f, 0.f);
    right = vec4(cross(direction.xyz(), up.xyz()).normalize(), 0.f);
    up = vec4(cross(right.xyz(), direction.xyz()).normalize(), 0.f);
    fov = 45;
    aspectRatio = 1;
    near = -0.5;
    far = -50;
    screen_w = 100;
    screen_h = 100;
}

Camera::Camera(vec3 pos, vec3 dir, vec3 u, double w, double h, double fov, double r, double n, double f) {
    position = vec4(pos, 1.f);
    direction = vec4(dir.normalize(), 0.f);
    right = vec4(cross(direction.xyz(), u).normalize(), 0.f);
    up = vec4(cross(right.xyz(), direction.xyz()).normalize(), 0.f);
    this->fov = fov*PI/180.0;
    aspectRatio = r;
    near = n;
    far = f;
    fov = fov;
    screen_w = w;
    screen_h = h;
}

void Camera::set(vec3 pos, vec3 dir, vec3 u, double w, double h, double fov, double r, double n, double f) {
    position = vec4(pos, 1.f);
    direction = vec4(dir.normalize(), 0.f);
    right = vec4(cross(direction.xyz(), u).normalize(), 0.f);
    up = vec4(cross(right.xyz(), direction.xyz()).normalize(), 0.f);
    this->fov = fov;
    aspectRatio = r;
    near = n;
    far = f;
    screen_w = w;
    screen_h = h;
}

double Camera::width() {
    return aspectRatio * height();
}

double Camera::height() {
    return -2*near*tan(fov*PI/180/2);
}

vec3 Camera::viewport_u() {
    return right.xyz().normalize() * width();
}

vec3 Camera::viewport_v() {
    return -up.xyz().normalize() * height();
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
    double h = -2*near*tan(fov*PI/180/2);
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