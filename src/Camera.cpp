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
    viewO = position.xyz();
    viewportO = viewO - direction.xyz() * near;
    viewportU = right.xyz() * width() / screen_w;
    viewportV = -up.xyz() * height() / screen_h;

    double defcusRadius = -near*tan(deg2rad(defocusAngle / 2));
    defocusDiskU = right.xyz() * defcusRadius;
    defocusDiskV = -up.xyz() * defcusRadius;
}

void Camera::render(const BVHNode &objects, TGAImage &image) {
    init();
    ProgressBar bar(10);
    vec3 color;

#pragma omp parallel for schedule(dynamic, 1) private(color)      // OpenMP
    for (int j = 0; j < screen_h; j ++) {
        std::string message = "Renderering (" + std::to_string(spp) + " spp)";
        bar.update(static_cast<double>(j)/(screen_h-1), message);
        std::vector<vec3> line;
        for (int i = 0; i < screen_w; i ++) {
            color = ZERO_VEC3;
            for (int k = 0; k < spp; k ++) {
                Ray ray = getRay(i, j);
                color += radiance(ray, 0, objects);
            }
            color = clampVec3(color / spp); // Clamp to [0, 1]
            gammaCorrection(color); // Gamma correction
            image.set(i, j, TGAColor(255*color));
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

inline vec2 Camera::pixSampleSquare() {
    // vec2 sample = sampleUnitSpuare() - 0.5;
    vec2 sample = 2.0*sampleUnitSpuare();
    // Tant filtering
    sample.x = sample.x < 1.0 ? sqrt(sample.x) - 1 : 1 - sqrt(sample.x);
    sample.y = sample.y < 1.0 ? sqrt(sample.y) - 1 : 1 - sqrt(sample.y);
    return sample;
}

inline vec2 Camera::pixSampleDisk(double radius) {
    vec2 sample = sampleUnitDisk();
    sample *= radius;
    return sample;
}

Ray Camera::getRay(const int &i, const int &j) {
    vec2 pixSamp = pixSampleSquare();
    vec3 pixPos = viewportO + (i - (screen_w >> 1) + 0.5)*viewportU + (j - (screen_h >> 1) + 0.5)*viewportV;
    pixPos += viewportU*pixSamp.x + viewportV*pixSamp.y;
    vec3 rayOrig = viewO;
    if (defocusAngle > 0) {
        vec2 defcusSamp = pixSampleDisk();
        rayOrig += defocusDiskU*defcusSamp.x + defocusDiskV*defcusSamp.y;
    }
    vec3 rayDir = (pixPos - rayOrig).normalize();
    return Ray(rayOrig, rayDir);
}

vec3 Camera::radiance(const Ray &r, int depth, const BVHNode &objects) {
    InterRecord rec;
    Interval rayt(1e-3, INF);

    if (!objects.intersect(r, rayt, rec))
        return backgroundLight;

    Ray scattered;
    vec3 attenuation;
    vec3 c = rec.mat->color();
    vec3 e = rec.mat->emit();
    double p = c.x>c.y && c.x>c.z ? c.x : c.y>c.z ? c.y : c.z; // max refl

    if (++depth > maxDepth) {
        // if (p > randDouble()) c /= p;
        // else return e;
        return ZERO_VEC3;
    }

    if (!rec.mat->scatter(r, rec, scattered, attenuation))
        return e;

    return e + attenuation.mult(radiance(scattered, depth, objects));
}
