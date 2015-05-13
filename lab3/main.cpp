#include "Lexer.h"
#include "Container.h"
#include "Simplex.h"

#include <iostream>
#include <fstream>
#include <cassert>
#include <algorithm>

#include <delaunay2.h>
#include "vtk.h"

using namespace config;

struct IntPoint {
    const Point *p;
    bool bnd;
    RhsVal rhs;
};

struct BndPoint {
    const Point *p;
    Point n;
    BcVal bc;
};

typedef Container<IntPoint> Interior;
typedef Container<BndPoint> Boundary;

std::ofstream logfile;

void mesh(const Region *r, const Problem *prob, Boundary &bnd, Interior &inter, double hbord, double h, double aspect, int lays) {
    const auto &cont = r->contour(hbord);

    for (const auto &p : cont) {
        if (!r->border(p))
            continue;

        const Point n = r->normal(p);
        const BcVal bc = r->condition(p);

        const Point *qp = new Point(p);
        bool added = bnd.add(BndPoint{qp, n, bc});
        added = inter.add(IntPoint{qp, true, prob->eval(p)});
        if (!added)
            delete qp;

        double s = aspect;

        for (int i = 1; i <= lays; i++) {
            const Point pp(p.x - s * i * hbord * n.x, p.y - s * i * hbord * n.y);
            if (r->inside(pp)) {
                const Point *qp = new Point(pp);
                bool added = inter.add(IntPoint{qp, false, prob->eval(pp)});
                if (!added)
                    delete qp;
            }
        }
    }

    const auto &bnds = r->bounds();
    Point ll = bnds.first;
    Point ur = bnds.second;

    ll.x -= h;
    ll.y -= h;
    ur.x -= h;
    ur.y -= h;

    int row = 0;
    for (double y = ll.y; y <= ur.y; y += 0.866025 * h) {
        for (double x = ll.x + 0.5 * row * h; x <= ur.x; x += h) {
            const Point p(x, y);
            if (r->inside(p)) {
                const Point *q = new Point(p);
                bool added = inter.add(IntPoint{q, false, prob->eval(p)});
                if (!added)
                    delete q;
            }
        }
        row = 1 - row;
    }

}

bool computeInternalEquation(
        const Problem *prob,
        const IntPoint &ipoint,
        const std::vector<IntPoint> &idata,
        std::vector<std::pair<int, double>> &neibs,
        std::vector<double> &alpha,
        double &alpha0, double &beta,
        std::vector<const Point *> &bad
)
{
    const Point *pc = ipoint.p;

    double e0  = prob->e0;
    double e1  = prob->e1;
    double e2  = prob->e2;
    double e11 = prob->e11;
    double e12 = prob->e12;
    double e22 = prob->e22;

    double f  = ipoint.rhs.f;
    double fx = ipoint.rhs.fx;
    double fy = ipoint.rhs.fy;

    const int K = 7;

    std::vector<double> A(K * neibs.size());
    std::vector<double> b(K);
    std::vector<double> c(neibs.size());

    int k = 0;
    for (const auto &jv : neibs) {
        int j = jv.first;
        const Point &p = *idata[j].p;
        double x = p.x - pc->x;
        double y = p.y - pc->y;
        A[K * k + 0] = x - e0 / e11 * x * x * x / 6;
        A[K * k + 1] = y - e0 / e22 * y * y * y / 6;
        A[K * k + 2] = x * x / 2 - e1 / e11 * x * x * x / 6;
        A[K * k + 3] = x * y - e2 / e11 * x * x * x / 6 - e1 / e22 * y * y * y / 6;
        A[K * k + 4] = y * y / 2 - e2 / e22 * y * y * y / 6;
        A[K * k + 5] = x * x * y / 2 - (e12 / e11 * x * x * x + e11 / e22 * y * y * y) / 6;
        A[K * k + 6] = y * y * x / 2 - (e22 / e11 * x * x * x + e12 / e22 * y * y * y) / 6;
        c[k] = pow(x * x + y * y, 1.5);
        k++;
    }
    b[0] = e1;
    b[1] = e2;
    b[2] = e11;
    b[3] = e12;
    b[4] = e22;
    b[5] = 0;
    b[6] = 0;

    std::vector<int> nonzero;
    std::vector<double> alp;
    simplex::ExtraStatus es;
    simplex::Status res = simplex::solveProblem(A, b, c, nonzero, alp, nullptr, &es);
    if (res != simplex::Status::Solved) {
        logfile << "Failed to solve problem at interior point " << *pc << ". Neibs = {" << std::endl;
        for (const auto &jv : neibs) {
            int j = jv.first;
            const Point &p = *idata[j].p;
            double x = p.x - pc->x;
            double y = p.y - pc->y;
            logfile << "{" << x << ", " << y << "}, " << std::endl;;
        }
        logfile << "}" << std::endl;
    }

    const auto oldneibs = neibs;
    neibs.clear();
    alpha.resize(neibs.size());

    if (res != simplex::Status::Solved) {
        if (res != simplex::Status::Infeasible || es != simplex::ExtraStatus::Normal) {
            logfile << "[I] Status = " << textStatus(res)
                << ", ExtraStatus = " << es.value() << std::endl;
        }
        bad.push_back(pc);
        return false;
    }

    for (const auto &v : nonzero)
        neibs.push_back(oldneibs[v]);

    alpha = alp;
    alpha0 = e0;

    beta = f;
    for (size_t j = 0; j < neibs.size(); j++) {
        const Point &p = *idata[neibs[j].first].p;
        double x = p.x - pc->x;
        double y = p.y - pc->y;

        alpha0 -= alp[j];
        beta += alp[j] * (x * x * x * fx / 6 / e11 + y * y * y * fy / 6 / e22);
    }

    return true;
}

bool computeBoundaryEquation(
        const Problem *prob, double f,
        const BndPoint &bdata,
        const std::vector<IntPoint> &idata,
        std::vector<std::pair<int, double>> &neibs,
        std::vector<double> &alpha,
        double &alpha0, double &beta,
        std::vector<const Point *> &bad
)
{
    const Point *pc = bdata.p;

    double nux = bdata.n.x;
    double nuy = bdata.n.y;

    double omega = bdata.bc.a;
    double sigma = bdata.bc.b;
    double rho = bdata.bc.g;

    double e0  = prob->e0;
    double e1  = prob->e1;
    double e2  = prob->e2;
    double e11 = prob->e11;
    double e12 = prob->e12;
    double e22 = prob->e22;

    double eps = 1e-6;

    if (omega <= -eps) {
        omega = -omega;
        sigma = -sigma;
        rho = -rho;
    }

    if (std::abs(omega) < eps && sigma < 0) {
        omega = 0;
        sigma = -sigma;
        rho = -rho;
    }

    const int K = 4;

    std::vector<double> A(K * neibs.size());
    std::vector<double> b(K);
    std::vector<double> c(neibs.size());

    int k = 0;
    for (const auto &jv : neibs) {
        int j = jv.first;
        const Point &p = *idata[j].p;
        double x = p.x - pc->x;
        double y = p.y - pc->y;
        A[K * k + 0] = x - 0.5 * e1 * (x * x + y * y) / (e11 + e22);
        A[K * k + 1] = y - 0.5 * e2 * (x * x + y * y) / (e11 + e22);
        A[K * k + 2] = x * y - 0.5 * e12 * (x * x + y * y) / (e11 + e22);
        A[K * k + 3] = 0.5 * (e22 * x * x - e11 * y * y);
        c[k] = pow(x * x + y * y, 1.5);
        k++;
    }

    b[0] = -omega * nux;
    b[1] = -omega * nuy;
    b[2] = 0;
    b[3] = 0;

    std::vector<int> nonzero;
    std::vector<double> alp;
    simplex::ExtraStatus es;
    simplex::Status res = simplex::solveProblem(A, b, c, nonzero, alp, nullptr, &es);
    if (res != simplex::Status::Solved) {
        logfile << "Failed to solve problem at boundary point " << *pc << std::endl;
        for (const auto &jv : neibs) {
            int j = jv.first;
            const Point &p = *idata[j].p;
            double x = p.x - pc->x;
            double y = p.y - pc->y;
            logfile << "Neib " << j << ": delta = " << Point(x, y) << std::endl;;
        }
    }

    const auto oldneibs = neibs;
    neibs.clear();
    alpha.resize(neibs.size());

    if (res != simplex::Status::Solved) {
        if (res != simplex::Status::Infeasible || es != simplex::ExtraStatus::Normal) {
            logfile << "[B] Status = " << textStatus(res)
                << ", ExtraStatus = " << es.value() << std::endl;
        }
        bad.push_back(pc);
        return false;
    }

    for (const auto &v : nonzero)
        neibs.push_back(oldneibs[v]);

    alpha = alp;
    alpha0 = -sigma;
    beta = -rho;

    for (size_t j = 0; j < neibs.size(); j++) {
        const Point &p = *idata[neibs[j].first].p;
        double x = p.x - pc->x;
        double y = p.y - pc->y;
        double gj = 0.5 / (e11 + e22) * (x * x + y * y);

        alpha0 -= alp[j];
        alpha0 += alp[j] * gj * e0;
        beta += alp[j] * gj * f;
    }

    return true;
}

std::vector<const Point *> buildCoeff(
    Interior inter, Boundary bnd, const Region *reg, const Problem *prob,
    std::vector<IntPoint> &idata,
    std::vector<int> &bndidx,
    std::vector<std::vector<std::pair<int, double> > > &neibs,
    std::vector<std::vector<double> > &alpha,
    std::vector<double> &alpha0,
    std::vector<double> &beta
)
{
    const auto bnds = reg->bounds();
    const double diam = bnds.second.x - bnds.first.x + bnds.second.y - bnds.first.y;

    const int nneib = 20;
    idata = inter.data();
    const auto &bdata = bnd.data();

    alpha.clear();
    alpha0.clear();
    beta.clear();
    neibs.clear();
    bndidx.clear();

    alpha.resize(idata.size());
    alpha0.resize(idata.size());
    beta.resize(idata.size());
    neibs.resize(idata.size());
    bndidx.assign(bdata.size(), -1);

    for (size_t i = 0; i < idata.size(); i++) {
        const Point pa(*idata[i].p);
        const auto &b = inter.lookupBucket(pa);

        for (size_t j = i + 1; j < idata.size(); j++) {
            const Point pb(*idata[j].p);
            if (!b.inall(pb))
                continue;
            double dd = b.distance(pa, pb);
            for (double s = 0.05; s < 1; s += 0.1) {
                Point mid(s * pa.x + (1 - s) * pb.x, s * pa.y + (1 - s) * pb.y);
                if (!reg->inside(mid))
                    dd = diam;
            }
            neibs[i].push_back(std::make_pair(j, dd));
            neibs[j].push_back(std::make_pair(i, dd));
        }

        bool isbnd = idata[i].bnd;
        if (!isbnd)
            continue;
        for (size_t j = 0; j < bdata.size(); j++)
            if (idata[i].p == bdata[j].p)
                bndidx[j] = i;
    }
    for (size_t i = 0; i < idata.size(); i++) {
        std::sort(neibs[i].begin(), neibs[i].end(),
            [](const std::pair<int, double> &a, const std::pair<int, double> &b) -> bool {
                return a.second < b.second;
            }
        );
        if (neibs[i].size() < nneib)
            std::cerr << "Point " << idata[i].p << " has only " << neibs[i].size() << " neibs!" << std::endl;
        else
            neibs[i].erase(neibs[i].begin() + nneib, neibs[i].end());
    }

    int badint = 0;
    int badbnd = 0;

    std::vector<const Point *> bad;

    for (size_t i = 0; i < neibs.size(); i++) {
        if (!idata[i].bnd) {
            bool f = computeInternalEquation(prob, idata[i], idata, neibs[i], alpha[i], alpha0[i], beta[i], bad);
            if (!f)
                badint++;
        }
    }

    for (size_t ib = 0; ib < bdata.size(); ib++) {
        size_t i = bndidx[ib];
        double rhsf = idata[i].rhs.f;
        bool f = computeBoundaryEquation(prob, rhsf, bdata[ib], idata, neibs[i], alpha[i], alpha0[i], beta[i], bad);
        if (!f)
            badbnd++;
    }

    std::cout << "Found " << badint << " bad internal points and " << badbnd << " bad boundary points" << std::endl;

    return bad;
}

int main(int argc, char **argv) {
	if (argc < 2) {
		std::cerr << "USAGE: ellipt <config.geom>" << std::endl;
		return 1;
	}

    logfile.open("bad_points.log");

	std::fstream f(argv[1], std::ios::in);

    if (!f) {
        std::cerr << "Could not open file `" << argv[1] << "'" << std::endl;
        return 1;
    }

	Lexer lexer(f, argv[1]);

	const Config *configraw = 0;
	Parser parser(lexer, configraw);
//    parser.set_debug_level(10);

	if (parser.parse())
		return 1;

    std::unique_ptr<const Config> config(configraw);

    const Region *reg = config->r;

    std::cout << reg->print() << std::endl;

    double hbord = config->ms.hbound;
    double h = config->ms.hinter;
    double aspect = config->ms.aspect;
    int lays = 3;

    const auto &bnds = reg->bounds();
    Point ll = bnds.first;
    Point ur = bnds.second;
    const double diam = ur.x - ll.x + ur.y - ll.y;

    ll.x -= 10 * h;
    ll.y -= 10 * h;
    ur.x += 10 * h;
    ur.y += 10 * h;

    const double q = 5;

    Interior inter(ll, ur, q * h);
    Boundary bnd(ll, ur, q * h);

    const Problem *prob = config->p;

    mesh(reg, prob, bnd, inter, hbord, h, aspect, lays);

    std::vector<IntPoint> idata;
    std::vector<int> bndidx;
    std::vector<std::vector<std::pair<int, double> > > neibs;
    std::vector<std::vector<double> > alpha;
    std::vector<double> alpha0;
    std::vector<double> beta;

    int maxiter = 10;

    for (int iter = 0; ; iter++) {
        const auto bad = buildCoeff(inter, bnd, reg, prob,
                idata, bndidx, neibs, alpha, alpha0, beta
            );
        std::cout << "From " << idata.size() << " points "<< bad.size() << " are marked as bad" << std::endl;
#if 0
        std::ofstream f("neibs." + std::to_string(iter));
        for (size_t i = 0; i < neibs.size(); i++)
            if (!idata[i].bnd)
                for (const auto &jv : neibs[i]) {
                    int j = jv.first;
                    f << idata[i].p->x << " " << idata[i].p->y << std::endl;
                    f << idata[j].p->x << " " << idata[j].p->y << std::endl << std::endl;
                }
        f.close();

        std::ofstream ff("bnds." + std::to_string(iter));
        for (size_t i = 0; i < neibs.size(); i++)
            if (idata[i].bnd)
                for (const auto &jv : neibs[i]) {
                    int j = jv.first;
                    ff << idata[i].p->x << " " << idata[i].p->y << std::endl;
                    ff << idata[j].p->x << " " << idata[j].p->y << std::endl << std::endl;
                }
        ff.close();

        std::ofstream g("bad." + std::to_string(iter));
        for (size_t i = 0; i < bad.size(); i++)
            g << bad[i]->x << " " << bad[i]->y << std::endl;
        g.close();
#endif
        if (bad.empty())
            break;

        for (const auto &p : bad) {
            inter.remove(p);
            bnd.remove(p);
        }

        if (iter >= maxiter) {
            std::cerr << "Could not remove all bad points in " << maxiter << " iterations. Try different stepsize" << std::endl;
            abort();
        }
    }

    ll.x -= diam / 2;
    ll.y -= diam / 2;
    ur.x += diam / 2;
    ur.y += diam / 2;

    std::vector<del_point2d_t> points;
    for (const auto &v : idata)
        points.push_back(del_point2d_t{v.p->x, v.p->y});

    // Add bounding box
    points.push_back(del_point2d_t{ll.x, ll.y});
    points.push_back(del_point2d_t{ur.x, ll.y});
    points.push_back(del_point2d_t{ur.x, ur.y});
    points.push_back(del_point2d_t{ll.x, ur.y});

    size_t numpts = idata.size();
#if 0
    double dmin = diam;
    const Bucket<IntPoint> &b = inter.bucks[0];
    int i0 = -1, j0 = -1;
    for (size_t i = 0; i < numpts; i++)
        for (size_t j = i + 1; j < numpts; j++) {
            double d = b.distance(Point(points[i].x, points[i].y), Point(points[j].x, points[j].y));
            if (d < dmin) {
                dmin = d;
                i0 = i;
                j0 = j;
            }
        }
    std::cout << "dmin = " << dmin << ", i0 = " << i0 << ", j0 = " << j0 << std::endl;

    std::ofstream fpts("points");
    fpts.precision(30);
    for (const auto &p : points) {
        fpts << p.x << " " << p.y << std::endl;
    }
    fpts.close();

    std::ifstream ipts("points");
    for (auto &p : points) {
        ipts >> p.x >> p.y;
    }
    ipts.close();
#endif

    delaunay2d_t *tri = delaunay2d_from(points.data(), points.size(), nullptr);
    std::vector<del_triface_t> trifilt;
    int k = 0;
    for (unsigned i = 0; i < tri->num_faces; i++) {
        unsigned nvert = tri->faces[k++];
        unsigned center = tri->faces[k++];
        unsigned next = tri->faces[k++];
        Point pc(points[center].x, points[center].y);
        for (unsigned j = 0; j < nvert - 2; j++) {
            unsigned prev = next;
            next = tri->faces[k++];
            Point pp(points[prev].x, points[prev].y);
            Point pn(points[next].x, points[next].y);
            Point c((pc.x + pp.x + pn.x) / 3, (pc.y + pp.y + pn.y) / 3);
            if (reg->inside(c) && center < numpts && prev < numpts && next < numpts)
                trifilt.push_back(del_triface_t{center, prev, next});
        }
    }
    delaunay2d_release(tri);

    points.pop_back();
    points.pop_back();
    points.pop_back();
    points.pop_back();

    /**
     *
     * uhat - u
     * -------- = alpha0 u + sum over neibs alpha uneib - beta
     *    dt
     *
     * 1 + alpha0 dt > 0
     * */

    double dt = 1e6;
    for (size_t i = 0; i < numpts; i++)
        if (1 + alpha0[i] * dt < 0)
            dt = -1 / alpha0[i];

    std::cout << "Maximum possible dt = " << dt << std::endl;

    std::vector<double> U(numpts, 0);
    std::vector<double> Unew(numpts, 0);

    for (const auto &v : idata)
        delete v.p;

    for (int iteration = 0; iteration < 10000000; iteration ++) {
        double diff = 0;
        for (size_t i = 0; i < numpts; i++) {
            Unew[i] = -dt * beta[i] + (1 + dt * alpha0[i]) * U[i];
            const auto &myneibs = neibs[i];
            for (size_t j = 0; j < myneibs.size(); j++)
                Unew[i] += dt * alpha[i][j] * U[myneibs[j].first];
            diff += std::abs(U[i] - Unew[i]);
        }

        std::swap(U, Unew);

        bool converged = diff < 1e-6;

        if ((iteration % 10000 == 0) || converged) {
            std::cout << "iteration " << iteration << ", difference = " << diff << std::endl;
            save(argv[1], iteration, points, trifilt, U);
        }

        if (converged)
            break;
    }

	return 0;
}
