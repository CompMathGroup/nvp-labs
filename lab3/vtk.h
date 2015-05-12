#ifndef __VTK_H__
#define __VTK_H__

#include <delaunay2.h>
#include <string>
#include <vector>

void save(const std::string &prefix, int step,
        const std::vector<del_point2d_t> &pts,
        const std::vector<del_triface_t> &faces,
        const std::vector<double> &U);

#endif
