#ifndef __TRIANGLE_H__
#define __TRIANGLE_H__

#include "Intersectable.h"

/******************************************************************
*                        Class Triangle                           *
******************************************************************/
struct Vertex {
    vec3 pos;
    vec3 norm;
    vec2 uv;

    Vertex(const vec3 &p) {
        pos = p;
    }
};

class Triangle: public Intersectable {
public:
    std::vector<Vertex> *verts;
    int vertIndex[3];
    shared_ptr<Material> mat;
    BoundingBox bbox;

    Triangle(std::vector<Vertex> *vertlist, int v1, int v2, int v3, shared_ptr<Material> mat=nullptr) {
        verts = vertlist;
        vertIndex[0] = v1;
        vertIndex[1] = v2;
        vertIndex[2] = v3;
        this->mat = mat;
        update();
    }

    Triangle(std::vector<Vertex> *vertlist, int index[3], shared_ptr<Material> mat=nullptr) {
        verts = vertlist;
        vertIndex[0] = index[0];
        vertIndex[1] = index[1];
        vertIndex[2] = index[2];
        this->mat = mat;
        update();
    }

    Vertex vertex(int index) const {
        return (*verts)[vertIndex[index]];
    }

    void update() override {
        bbox = boundingBox();
        u = vertex(1).pos - vertex(0).pos;
        v = vertex(2).pos - vertex(0).pos;
        n = u.cross(v);
        w = n / (n*n);
        D = n*vertex(0).pos;
    }

    void setMaterial(shared_ptr<Material> mat) override {
        this->mat = mat;
    }

    BoundingBox boundingBox() const override {
        auto ax = fmin(vertex(0).pos.x, fmin(vertex(1).pos.x, vertex(2).pos.x));
        auto ay = fmin(vertex(0).pos.y, fmin(vertex(1).pos.y, vertex(2).pos.y));
        auto az = fmin(vertex(0).pos.z, fmin(vertex(1).pos.z, vertex(2).pos.z));
        auto bx = fmax(vertex(0).pos.x, fmax(vertex(1).pos.x, vertex(2).pos.x));
        auto by = fmax(vertex(0).pos.y, fmax(vertex(1).pos.y, vertex(2).pos.y));
        auto bz = fmax(vertex(0).pos.z, fmax(vertex(1).pos.z, vertex(2).pos.z));
        return BoundingBox(vec3(ax, ay, az), vec3(bx, by, bz));
    }

    vec2 uv(const vec3 &p) {
        auto alpha = w*p.cross(v);
        auto beta = w*u.cross(p);
        return vec2(alpha, beta);
    }

    // Compute wheather a ray hitted the Triangle
    bool intersect(const Ray &ray, Interval rayt, InterRecord &rec) const override {
        double k = n*ray.dir;
        // if (k > -EPS) return false;    // No intersection if triangle is back-face
        if (fabs(k) < EPS) return false;    // No intersection if triangle is back-face
        double t = (D - ray.orig*n) / k;
        if (!rayt.contain(t)) return false;
        vec3 intersection = ray.at(t);
        vec3 p = intersection - vertex(0).pos;
        auto alpha = w*p.cross(v);
        auto beta = w*u.cross(p);
        if (alpha < 0 || beta < 0 || (1 - alpha - beta) < 0) return false;

        // Update rec
        vec3 normal = n.normalize();
        rec.t = t;
        rec.p = intersection;
        rec.uv = vec2(alpha, beta);
        rec.inward = (k < 0);
        rec.normal = rec.inward ? normal : -normal; 
        rec.mat = mat;

        return true;
    }

private:
    double D;
    vec3 u, v, w, n;

};

#endif