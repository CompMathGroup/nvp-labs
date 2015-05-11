#ifndef __CONFIG__VALUE_H__
#define __CONFIG__VALUE_H__

#include <string>
#include <vector>
#include <ostream>
#include <cmath>
#include <cassert>

namespace config {

struct Point {
    double x, y;
    Point(double x, double y) : x(x), y(y) { }
};

inline std::ostream &operator<<(std::ostream &o, const Point &p) {
    return o << "(" << p.x << ", " << p.y << ")";
}

/* a du/dn + b u = g */
struct BcVal {
    double a, b, g;
};

struct BoundaryCondition {
    const std::string a, b, g;
    BoundaryCondition(const std::string &a, const std::string &b, const std::string &g)
        : a(a), b(b), g(g)
    {
        compile();
    }
    void compile();
    BcVal eval(const Point &p) const;
    virtual ~BoundaryCondition() { }
};

struct Diriclet : public BoundaryCondition {
    Diriclet(const std::string &g)
        : BoundaryCondition("return 0;", "return 1;", g)
    { }
};

struct Robin : public BoundaryCondition {
    Robin(const std::string &b, const std::string &g)
        : BoundaryCondition("return 1;", b, g)
    { }
};

struct Neumann : public Robin {
    Neumann(const std::string &g)
        : Robin("return 0;", g)
    { }
};

struct Path {
    std::vector<std::pair<const Point, const BoundaryCondition *> > path;
    Path() { }
    void append(const Point &p, const BoundaryCondition *bc) {
        path.push_back(std::make_pair(p, bc));
    }
};

struct Region {
    static constexpr double tol = 1e-6;
    virtual bool inside(const Point &p) const = 0;
    virtual bool outside(const Point &p) const = 0;
    virtual bool border(const Point &p) const {
        return !inside(p) && !outside(p);
    }
    virtual Point normal(const Point &p) const = 0;
    virtual BcVal condition(const Point &p) const = 0;
    virtual std::vector<Point> contour(double h) const = 0;
    virtual std::string print() const = 0;
    virtual std::pair<Point, Point> bounds() const = 0;
    virtual ~Region() { }
};

struct Polygon : public Region {
    std::vector<const BoundaryCondition *> bc;
    std::vector<Point> ps;
    std::vector<Point> psinner;
    std::vector<Point> psouter;

    Polygon(const Path &path);
    virtual bool inside(const Point &p) const override;
    virtual bool outside(const Point &p) const override;
    virtual Point normal(const Point &p) const override;
    virtual BcVal condition(const Point &p) const override;
    virtual std::vector<Point> contour(double h) const override;
    virtual std::pair<Point, Point> bounds() const override;
    virtual std::string print() const override {
        std::string ret = "Poly[";
        for (const auto &p : ps) {
            ret += "{";
            ret += std::to_string(p.x);
            ret += ", ";
            ret += std::to_string(p.y);
            ret += "}, ";
        }
        ret[ret.size() - 1] = ']';
        return ret;
    }
};

struct Circle : public Region {
    const Point center;
    const double radius;
    const BoundaryCondition *bc;
    Circle(const Point &p, const double r, const BoundaryCondition *bc)
        : center(p), radius(r)
    { }
    virtual bool inside(const Point &p) const override {
        double dx = p.x - center.x;
        double dy = p.y - center.y;
        return dx * dx + dy * dy <= radius * radius * (1 - 2 * tol);
    }
    virtual bool outside(const Point &p) const override {
        double dx = p.x - center.x;
        double dy = p.y - center.y;
        return dx * dx + dy * dy >= radius * radius * (1 + 2 * tol);
    }
    virtual Point normal(const Point &p) const override {
        double dx = p.x - center.x;
        double dy = p.y - center.y;
        double dr = std::sqrt(dx * dx + dy * dy);
        return Point(dx / dr, dy / dr);
    }
    virtual BcVal condition(const Point &p) const override {
        return bc->eval(p);
    }
    virtual std::string print() const override {
        return "Circle[{" + std::to_string(center.x) + ", " +
            std::to_string(center.y) + "}, " + std::to_string(radius) + "]";
    }
    virtual std::vector<Point> contour(double h) const override {
        const double pi = 4 * std::atan(1.);
        int N = 2 * pi * radius / h;
        double dphi = 2 * pi / N;
        if (N < 6)
            N = 6;
        std::vector<Point> ret;
        for (int i = 0; i < N; i++) {
            ret.push_back(Point(
                center.x + radius * cos(i * dphi),
                center.y + radius * sin(i * dphi)
            ));
        }
        return ret;
    }
    virtual std::pair<Point, Point> bounds() const override {
        return std::make_pair(
            Point(center.x - radius, center.y - radius),
            Point(center.x + radius, center.y + radius)
        );
    }
};

struct Union : public Region {
    const Region *r1, *r2;
    Union(const Region *r1, const Region *r2) : r1(r1), r2(r2) { }
    virtual bool inside(const Point &p) const override {
        return r1->inside(p) || r2->inside(p);
    }
    virtual bool outside(const Point &p) const override {
        return r1->outside(p) && r2->outside(p);
    }
    virtual Point normal(const Point &p) const override {
        assert(border(p));
        if (r1->border(p))
            return r1->normal(p);
        return r2->normal(p);
    }
    virtual BcVal condition(const Point &p) const override {
        assert(border(p));
        if (r1->border(p))
            return r1->condition(p);
        return r2->condition(p);
    }
    virtual std::string print() const override {
        return "(" + r1->print() + " + " + r2->print() + ")";
    }
    virtual std::vector<Point> contour(double h) const override {
        std::vector<Point> all(r1->contour(h));
        const auto &add = r2->contour(h);
        all.insert(all.end(), add.begin(), add.end());
        return all;
    }
    virtual std::pair<Point, Point> bounds() const override {
        const auto &b1 = r1->bounds();
        const auto &b2 = r2->bounds();
        return std::make_pair(
            Point(std::min(b1.first.x, b2.first.x), std::min(b1.first.y, b2.first.y)),
            Point(std::max(b1.second.x, b2.second.x), std::max(b1.second.y, b2.second.y))
        );
    }
};

struct Subtract : public Region {
    const Region *r1, *r2;
    Subtract(const Region *r1, const Region *r2) : r1(r1), r2(r2) { }
    virtual bool inside(const Point &p) const override {
        return r1->inside(p) && r2->outside(p);
    }
    virtual bool outside(const Point &p) const override {
        return r1->outside(p) || r2->inside(p);
    }
    virtual Point normal(const Point &p) const override {
        assert(border(p));
        if (r1->border(p))
            return r1->normal(p);
        Point n(r2->normal(p));
        return Point(-n.x, -n.y);
    }
    virtual BcVal condition(const Point &p) const override {
        assert(border(p));
        if (r1->border(p))
            return r1->condition(p);
        return r2->condition(p);
    }
    virtual std::string print() const override {
        return "(" + r1->print() + " - " + r2->print() + ")";
    }
    virtual std::vector<Point> contour(double h) const override {
        std::vector<Point> all(r1->contour(h));
        const auto &add = r2->contour(h);
        all.insert(all.end(), add.begin(), add.end());
        return all;
    }
    virtual std::pair<Point, Point> bounds() const override {
        const auto &b1 = r1->bounds();
        const auto &b2 = r2->bounds();
        return std::make_pair(
            Point(std::min(b1.first.x, b2.first.x), std::min(b1.first.y, b2.first.y)),
            Point(std::max(b1.second.x, b2.second.x), std::max(b1.second.y, b2.second.y))
        );
    }
};

struct Intersect : public Region {
    const Region *r1, *r2;
    Intersect(const Region *r1, const Region *r2) : r1(r1), r2(r2) { }
    virtual bool inside(const Point &p) const override {
        return r1->inside(p) && r2->inside(p);
    }
    virtual bool outside(const Point &p) const override {
        return r1->outside(p) || r2->outside(p);
    }
    virtual Point normal(const Point &p) const override {
        assert(border(p));
        if (r1->border(p))
            return r1->normal(p);
        return r2->normal(p);
    }
    virtual BcVal condition(const Point &p) const override {
        assert(border(p));
        if (r1->border(p))
            return r1->condition(p);
        return r2->condition(p);
    }
    virtual std::string print() const override {
        return "(" + r1->print() + " * " + r2->print() + ")";
    }
    virtual std::vector<Point> contour(double h) const override {
        std::vector<Point> all(r1->contour(h));
        const auto &add = r2->contour(h);
        all.insert(all.end(), add.begin(), add.end());
        return all;
    }
    virtual std::pair<Point, Point> bounds() const override {
        const auto &b1 = r1->bounds();
        const auto &b2 = r2->bounds();
        return std::make_pair(
            Point(std::min(b1.first.x, b2.first.x), std::min(b1.first.y, b2.first.y)),
            Point(std::max(b1.second.x, b2.second.x), std::max(b1.second.y, b2.second.y))
        );
    }
};

}

#endif
