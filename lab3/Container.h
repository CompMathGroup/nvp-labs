#ifndef __CONTAINER_H__
#define __CONTAINER_H__

#include "Value.h"

#include <list>

namespace config {

template<class T>
struct Bucket {
    Point c;
    double a;
    std::list<T> data;

    Bucket() : c(0, 0), a(0) { }
    bool incore(const Point &p) {
        bool xin = (p.x > c.x - 0.5 * a) && (p.x <= c.x + 0.5 * a);
        bool yin = (p.y > c.y - 0.5 * a) && (p.y <= c.y + 0.5 * a);
        return xin && yin;
    }
    bool inall(const Point &p) {
        bool xin = (p.x > c.x - a) && (p.x <= c.x + a);
        bool yin = (p.y > c.y - a) && (p.y <= c.y + a);
        return xin && yin;
    }
    double distance(const Point &p1, const Point &p2) {
        return std::abs(p2.x - p1.x) + std::abs(p2.y - p1.y);
    }
    typename std::list<T>::iterator find_closest(const Point &p, double &dist) {
        dist = distance(p, data.front().p);
        typename std::list<T>::iterator ret = data.begin();
        for (auto it = data.begin(); it != data.end(); ++it) {
            double dd = distance(p, it->p);
            if (dd < dist) {
                dist = dd;
                ret = it;
            }
        }
        return ret;
    }
    void add(const T &v) {
        data.push_back(v);
    }
    void remove(const Point &p, double da) {
        if (!inall(p))
            return;
        data.remove_if([this, &p, da](const T &v) -> bool { return this->distance(p, v.p) < da; });
    }
};

template<class T>
struct Container {
    std::vector<Bucket<T> > bucks;
    double merge_tol;

    Container(const Point &ll, const Point &ur, const double bucketSize, const double merge_tol = 1e-3)
        : merge_tol(merge_tol)
    {
        double dx = ur.x - ll.x;
        double dy = ur.y - ll.y;
        const int nx = 2 + dx / bucketSize;
        const int ny = 2 + dy / bucketSize;
        bucks.resize(nx * ny);

        double a = bucketSize;

        for (int i = 0; i < nx; i++)
            for (int j = 0; j < ny; j++) {
                Bucket<T> b;
                b.a = bucketSize;
                b.c = Point(ll.x + i * a, ll.y + j * a);
                bucks.push_back(b);
            }
    }
    void add(const T &v) {
        for (auto &b : bucks)
            if (b.incore(v.p)) {
                double dist;
                const auto &w = b.find_closest(v.p, dist);
                if (dist < b.a * merge_tol) {
                    std::cout << "Point " << v.p << " was merged with " << w->p << std::endl;
                    return;
                }
            }
        for (auto &b : bucks)
            if (b.inall(v.p))
                b.add(v);
    }
    void add_nomerge(const T &v) {
        for (auto &b : bucks)
            if (b.inall(v.p))
                b.add(v);
    }
    void remove(const Point &p) {
        for (auto &b : bucks)
            b.remove(p, b.a * merge_tol);
    }
    std::vector<std::pair<T, Bucket<T> *> > data() {
        std::vector<std::pair<T, Bucket<T> *> > ret;
        for (auto &b : bucks)
            for (const auto &v : b.data)
                if (b.incore(v.p))
                    ret.push_back(std::make_pair(v, &b));
        return ret;
    }
};

}

#endif
