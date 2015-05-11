#include "Value.h"

#include <iostream>

namespace config {

Polygon::Polygon(const Path &path) {
    for (const auto &it : path.path) {
        ps.push_back(it.first);
        bc.push_back(it.second);
    }
    ps.push_back(path.path[0].first);
    ps.push_back(path.path[1].first);
    for (size_t i = 1; i < ps.size() - 1; i++) {
        const Point &p =  ps[i];
        const Point &pp = ps[i-1];
        const Point &pn = ps[i+1];
        double d1x = p.x - pp.x;
        double d1y = p.y - pp.y;
        double d2x = pn.x - p.x;
        double d2y = pn.y - p.y;

        double n1x = d1y * tol;
        double n1y = -d1x * tol;
        double n2x = d2y * tol;
        double n2y = -d2x * tol;

        double n1s = n1x * n1x + n1y * n1y;
        double n2s = n2x * n2x + n2y * n2y;
        double n1n2 = n1x * n2x + n1y * n2y;
        double det = n1s * n2s - n1n2 * n1n2;
        double alpha = (n1s * n2s - n1n2 * n2s) / det;
        double beta = (n1s * n2s - n1n2 * n1s) / det;

        psouter.push_back(Point(
            p.x + alpha * n1x + beta * n2x,
            p.y + alpha * n1y + beta * n2y
        ));
        psinner.push_back(Point(
            p.x - alpha * n1x - beta * n2x,
            p.y - alpha * n1y - beta * n2y
        ));
    }
    psinner.push_back(psinner[0]);
    psouter.push_back(psouter[0]);
}

std::pair<Point, Point> Polygon::bounds() const {
    Point ll(ps[0]), ur(ps[0]);
    for (const auto &p : ps) {
        if (ll.x > p.x)
            ll.x = p.x;
        if (ur.x < p.x)
            ur.x = p.x;
        if (ll.y > p.y)
            ll.y = p.y;
        if (ur.y < p.y)
            ur.y = p.y;
    }
    return std::make_pair(ll, ur);
}

static inline double isLeft(const Point &P0, const Point &P1, const Point &P2) {
    return ((P1.x - P0.x) * (P2.y - P0.y)
          - (P2.x - P0.x) * (P1.y - P0.y));
}

static bool insideTest(const Point &p, const std::vector<Point> &ps) {
    int wn = 0;

    for (size_t i = 1; i < ps.size(); i++) {
        if (ps[i-1].y <= p.y) {
            if (ps[i].y > p.y)
                if (isLeft(ps[i-1], ps[i], p) > 0)
                    ++wn;
        } else {
            if (ps[i].y <= p.y)
                if (isLeft(ps[i-1], ps[i], p) < 0)
                    --wn;
        }
    }
    return wn != 0;
}

bool Polygon::inside(const Point &p) const {
    return insideTest(p, psinner);
}

bool Polygon::outside(const Point &p) const {
    return !insideTest(p, psouter);
}

Point Polygon::normal(const Point &p) const {
    assert(border(p));
    for (size_t i = 0; i < bc.size(); i++) {
        double dx = ps[i+1].x - ps[i].x;
        double dy = ps[i+1].y - ps[i].y;
        double dr2 = dx * dx + dy * dy;
        if (isLeft(ps[i], ps[i+1], p) < tol * dr2) {
            double dr = std::sqrt(dr2);
            return Point(dy / dr, -dx / dr);
        }
    }
    return Point(0, 0);
}

BcVal Polygon::condition(const Point &p) const {
    assert(border(p));
    for (size_t i = 0; i < bc.size(); i++) {
        double dx = ps[i+1].x - ps[i].x;
        double dy = ps[i+1].y - ps[i].y;
        double dr2 = dx * dx + dy * dy;
        if (isLeft(ps[i], ps[i+1], p) < tol * dr2)
            return bc[i]->eval(p);
    }
    return BcVal();
}

std::vector<Point> Polygon::contour(double h) const {
    std::vector<Point> ret;
    for (size_t i = 0; i < bc.size(); i++) {
        double ax = ps[i  ].x;
        double ay = ps[i  ].y;
        double bx = ps[i+1].x;
        double by = ps[i+1].y;
        double dx = bx - ax;
        double dy = by - ay;
        double dr = sqrt(dx * dx + dy * dy);
        int N = dr / h;
        if (N < 3)
            N = 3;
        for (int j = 0; j < N; j++) {
            double s = 1. * j / N;
            ret.push_back(Point(ax + dx * s, ay + dy * s));
        }
    }
    return ret;
}

}
