#include "Value.h"

#include <iostream>

namespace config {

void BoundaryCondition::compile() {
/*    std::cout << "a = " << a << std::endl;
    std::cout << "b = " << b << std::endl;
    std::cout << "g = " << g << std::endl; */
}

BcVal BoundaryCondition::eval(const Point &p) const {
 /*   std::cout << "eval at (" << p.x << ", " << p.y << ")" << std::endl; */
    return BcVal{1, 2, 3};
}

void Problem::compile() {
}

RhsVal Problem::eval(const Point &p) const {
    return RhsVal{1, 0, 0};
}

}
