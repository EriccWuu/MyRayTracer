#ifndef UTILITY_H
#define UTILITY_H

#include <memory>
#include <vector>
#include "Intersectable.h"
#include "Ray.h"
#include "Material.h"

// Usings
using std::shared_ptr;
using std::make_shared;
using std::vector;

typedef std::shared_ptr<Intersectable> IntersectablePtr;
typedef std::vector<IntersectablePtr> Interlist;
typedef std::shared_ptr<Material> Matptr;

#endif