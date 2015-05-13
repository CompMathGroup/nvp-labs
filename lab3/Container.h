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

    bool inall(const Point &p) const {
        bool xin = (p.x > c.x - a) && (p.x <= c.x + a);
        bool yin = (p.y > c.y - a) && (p.y <= c.y + a);
        return xin && yin;
    }
    double distance(const Point &p1, const Point &p2) const {
        return std::abs(p2.x - p1.x) + std::abs(p2.y - p1.y);
    }
    typename std::list<T>::const_iterator find_closest(const Point &p, double &dist) const {
        dist = 10 * a;
        typename std::list<T>::const_iterator ret = data.begin();
        for (auto it = data.begin(); it != data.end(); ++it) {
            double dd = distance(p, *it->p);
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
    void remove(const Point *p) {
        data.remove_if([this, p](const T &v) -> bool { return p == v.p; });
    }
};

template<class T>
struct Container {
    const Point ll, ur;
    std::vector<Bucket<T> > bucks;
    double merge_tol;
    int nx, ny;

    Container(const Point &ll, const Point &ur, const double bucketSize, const double merge_tol = 1e-3)
        : ll(ll), ur(ur), merge_tol(merge_tol)
    {
        double dx = ur.x - ll.x;
        double dy = ur.y - ll.y;
        nx = 2 + dx / bucketSize;
        ny = 2 + dy / bucketSize;

        double a = bucketSize;

        for (int i = 0; i < nx; i++)
            for (int j = 0; j < ny; j++)
                bucks.push_back(Bucket<T>{Point(ll.x + i * a, ll.y + j * a), bucketSize});

//        std::cout << "Created " << nx << " x " << ny << " buckets" << std::endl;
    }
    size_t lookupBucketIdx(const Point &p) const {
        double a = bucks[0].a;
//        std::cout << p << " " << ll << " " << a << std::endl;
        double di = (p.x - ll.x) / a;
        double dj = (p.y - ll.y) / a;
        di += .5;
        dj += .5;
        int i = di;
        int j = dj;
        size_t ret = i * ny + j;
//        std::cout << "For point " << p << " bucket id = (" << i << ", " << j << ") = " << ret << std::endl;
        return ret;
    }
    Bucket<T> &lookupBucket(const Point &p) {
        return bucks[lookupBucketIdx(p)];
    }
    const Bucket<T> &lookupBucket(const Point &p) const {
        return bucks[lookupBucketIdx(p)];
    }
    void add(const T &v) {
        auto &b = lookupBucket(*v.p);
        double dist;
        const auto &w = b.find_closest(*v.p, dist);
        if (dist < b.a * merge_tol) {
            std::cout << "Point " << *v.p << " was merged with " << *w->p << std::endl;
            return;
        }
        for (auto &b : bucks)
            if (b.inall(*v.p))
                b.add(v);
    }
    void remove(const Point *p) {
        for (auto &b : bucks)
            b.remove(p);
    }
    std::vector<T> data() const {
        std::vector<T> ret;
        for (size_t ib = 0; ib < bucks.size(); ib++) {
            const auto &b = bucks[ib];
            for (const auto &v : b.data) {
                if (ib == lookupBucketIdx(*v.p))
                    ret.push_back(v);
            }
        }

        return ret;
    }
};

}

#endif
