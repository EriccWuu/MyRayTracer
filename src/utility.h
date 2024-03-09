#ifndef UTILITY_H
#define UTILITY_H

#include "stdio.h"
#include "string"
#include <memory>
#include <vector>
#include "Intersectable.h"
#include "BVH.h"
#include "Ray.h"
#include "Material.h"

// Usings
using std::shared_ptr;
using std::make_shared;
using std::vector;

// typedef std::vector<std::shared_ptr<Intersectable>> Interlist;
typedef std::shared_ptr<Material> Matptr;

class ProgressBar {
public:
    explicit ProgressBar(int width) : width_(width) {}

    void update(double progress, const std::string &message) {
        int barWidth = static_cast<int>(progress * width_);
        std::string barstr = "[";
        for (int i = 0; i < width_; ++i)
            barstr += (i < barWidth) ? '=' : ' ';
        barstr += "]";
        fprintf(stderr, "\r%s: %s %6.2f%% ",message.c_str(), barstr.c_str(), progress*100.0);
        fflush(stderr);
    }

private:
    int width_;
};

#endif