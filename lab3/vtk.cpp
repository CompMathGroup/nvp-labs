#include "vtk.h"

#include <fstream>
#include <iostream>
#include <algorithm>

template<typename T>
static void put(std::ofstream &f, const T value) {
    union {
        char buf[sizeof(T)];
        T val;
    } helper;
    helper.val = value;
    std::reverse(helper.buf, helper.buf + sizeof(T));
    f.write(helper.buf, sizeof(T));
}

void save(const std::string &prefix, int step,
        const std::vector<del_point2d_t> &pts,
        const std::vector<del_triface_t> &faces,
        const std::vector<double> &U)
{
    std::string filename(prefix + "." + std::to_string(step) + ".vtk");
    std::ofstream f(filename, std::ios::binary);

    if (!f) {
        std::cerr << "Could not write output file `" + filename+ "'" << std::endl;
        return;
    }

    f << "# vtk DataFile Version 3.0\nOutput\nBINARY\nDATASET UNSTRUCTURED_GRID";
    f << "\nPOINTS " << pts.size() << " float\n";
    for (const auto &p : pts) {
        put<float>(f, p.x);
        put<float>(f, p.y);
        put<float>(f, 1);
    }
    f << "\nCELLS " << faces.size() << " " << 4 * faces.size() << "\n";
    for (const auto &t : faces) {
        unsigned int cnt = 3;
        put(f, cnt);
        put(f, t.i);
        put(f, t.j);
        put(f, t.k);
    }
    f << "\nCELL_TYPES " << faces.size() << "\n";
    const int VTK_TRIANGLE = 5;
    for (size_t i = 0; i < faces.size(); i++)
        put(f, VTK_TRIANGLE);
    f << "\nPOINT_DATA " << U.size();
    f << "\nSCALARS U float\nLOOKUP_TABLE default\n";
    for (const auto &u : U)
        put<float>(f, u);

    f.close();
}
