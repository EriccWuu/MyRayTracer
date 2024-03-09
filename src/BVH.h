#ifndef __BVH_H__
#define __BVH_H__

#include "Intersectable.h"

class BVHNode: public Intersectable {
public:
    BVHNode(const std::vector<shared_ptr<Intersectable>> &ordered_obj, size_t start, size_t end) {
        // [start, end)
        size_t n =  end - start;
        if (n == 1) {
            left = right = ordered_obj[start];
            // std::cout << start << '\n';
        }
        else if (n == 2) {
            left = ordered_obj[start];
            right = ordered_obj[start+1];
            // std::cout << start << ' ' << start + 1 << '\n';
        }
        else {
            size_t mid = start + n / 2;
            left = make_shared<BVHNode>(ordered_obj, start, mid);   // [start, mid)
            right = make_shared<BVHNode>(ordered_obj, mid, end);    // [mid, end)
        }

        aabb = BoundingBox(left->boundingBox(), right->boundingBox());
    }

    bool intersect(const Ray &ray, Interval rayt, InterRecord & rec) const override {
        if (!aabb.intersect(ray))
            return false;
        bool hit_left = left->intersect(ray, rayt, rec);
        bool hit_right = right->intersect(ray, Interval(rayt.min, hit_left ? rec.t : rayt.max), rec);

        return hit_left || hit_right;
    }

    BoundingBox boundingBox() const override { return aabb; }

private:
    shared_ptr<Intersectable> left;
    shared_ptr<Intersectable> right;
    BoundingBox aabb;
};


#endif