#ifndef __MESH_H__
#define __MESH_H__

#include "Intersectable.h"
#include "Triangle.h"
#include "BVH.h"

/******************************************************************
*                        Class Mesh                               *
******************************************************************/
class Mesh: public Intersectable {
public:
    shared_ptr<Material> mat;

    Mesh(std::vector<Vertex> &vertlist, std::vector<int> &index, shared_ptr<Material> mat=nullptr) {
        verts = vertlist;
        vertIndex = index;
        this->mat = mat;
        if (index.size() % 3 != 0) std::cerr << "Mesh ERROR: Vertices number error.\n";
        int ind[3];
        for (int i = 0; i < index.size(); i += 3) {
            ind[0] = index[i];
            ind[1] = index[i+1];
            ind[2] = index[i+2];
            triangleList.push_back(std::make_shared<Triangle>(&verts, ind, mat));
        }
        // Build BVH of Mesh object
        std::sort(triangleList.begin(), triangleList.end(), [&](const shared_ptr<Intersectable> &a, const shared_ptr<Intersectable> &b) {return compare(a, b);});
        bvh = BVHNode(triangleList, 0, triangleList.size());
    }

    void setMaterial(shared_ptr<Material> mat) override {
        for (auto triangle: triangleList)
            triangle->setMaterial(mat);
    }

    void transform(const mat4 &M) {
        for (auto &vert: verts) 
            vert.pos = (M*vec4(vert.pos, 1.0)).xyz();
        for (auto triangle: triangleList)
            triangle->update();
        // Build BVH of Mesh object
        std::sort(triangleList.begin(), triangleList.end(), [&](const shared_ptr<Intersectable> &a, const shared_ptr<Intersectable> &b) {return compare(a, b);});
        bvh = BVHNode(triangleList, 0, triangleList.size());
    }

    BoundingBox boundingBox() const override {
        return bvh.boundingBox();
    }

    // Compute wheather a ray hitted the Triangle
    bool intersect(const Ray &ray, Interval rayt, InterRecord &rec) const override {
        if (!bvh.intersect(ray, rayt, rec)) return false;
        return true;
    }

private:
    std::vector<Vertex> verts;
    std::vector<int> vertIndex;
    std::vector<shared_ptr<Intersectable>> triangleList;
    BVHNode bvh;

};

static Mesh rectangle(double a, double h) {
    vec3 Q1 = {-0.5*a, 0, -0.5*a}, Q2 = {0.5*a, h, 0.5*a};
    std::vector<Vertex> vertices = {Vertex(Q1), Vertex(Q1 + a*XA), Vertex(Q1 + a*ZA), Vertex(Q1 + a*(XA + ZA))
                                , Vertex(Q2), Vertex(Q2 - a*XA), Vertex(Q2 - a*ZA), Vertex(Q2 - a*(XA + ZA))};
    std::vector<int> index = {0, 1, 2, 1, 3, 2  // buttom
                            , 0, 2, 7, 2, 5, 7  // left
                            , 0, 7, 1, 7, 6, 1  // back
                            , 4, 6, 5, 6, 7, 5  // up
                            , 4, 3, 6, 3, 1, 6  // right 
                            , 4, 5, 2, 2, 3, 4  // front 
                            };
    return Mesh(vertices, index);
}

#endif