#include "Lexer.h"
#include "Container.h"
#include "Simplex.h"

#include <iostream>
#include <fstream>
#include <cassert>
#include <algorithm>

using namespace config;

struct IntPoint {
    Point p;
    bool bnd;
};

struct BndPoint {
    Point p, n;
    BcVal bc;
};

typedef Container<IntPoint> Interior;
typedef Container<BndPoint> Boundary;

std::ofstream logfile;

void mesh(const Region *r, Boundary &bnd, Interior &inter, double hbord, double h, int lays) {
    const auto &cont = r->contour(hbord);

    for (const auto &p : cont) {
        if (!r->border(p))
            continue;

        const Point n = r->normal(p);
        const BcVal bc = r->condition(p);

        bnd.add(BndPoint{p, n, bc});
        inter.add_nomerge(IntPoint{p,true});

        double s = 1;

        for (int i = 1; i <= lays; i++) {
            const Point pp(p.x - s * i * hbord * n.x, p.y - s * i * hbord * n.y);
            if (r->inside(pp))
                inter.add_nomerge(IntPoint{pp,false});
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
            Point p(x, y);
            if (r->inside(p))
                inter.add_nomerge(IntPoint{p,false});
        }
        row = 1 - row;
    }

}

bool computeInternalEquation(
        const Point &pc,
        const std::vector<std::pair<IntPoint, Bucket<IntPoint> *>> &idata,
        std::vector<std::pair<int, double>> &neibs,
        std::vector<double> &alpha,
        double &alpha0, double &beta,
        std::vector<Point> &bad
)
{
    double e0 = 0;
    double e1 = 0;
    double e2 = 0;
    double e11 = 1;
    double e12 = 0;
    double e22 = 1;

    double f = 1;
    double fx = 0;
    double fy = 0;

    const int K = 7;

    std::vector<double> A(K * neibs.size());
    std::vector<double> b(K);
    std::vector<double> c(neibs.size());

    int k = 0;
    for (const auto &jv : neibs) {
        int j = jv.first;
        const Point &p = idata[j].first.p;
        double x = p.x - pc.x;
        double y = p.y - pc.y;
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
        logfile << "Failed to solve problem at interior point " << pc << ". Neibs = {" << std::endl;
        for (const auto &jv : neibs) {
            int j = jv.first;
            const Point &p = idata[j].first.p;
            double x = p.x - pc.x;
            double y = p.y - pc.y;
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
        const Point &p = idata[neibs[j].first].first.p;
        double x = p.x - pc.x;
        double y = p.y - pc.y;

        alpha0 -= alp[j];
        beta += alp[j] * (x * x * x * fx / 6 / e11 + y * y * y * fy / 6 / e22);
    }

    return true;
}

bool computeBoundaryEquation(
        const Point &pc,
        const std::vector<std::pair<IntPoint, Bucket<IntPoint> *>> &idata,
        const BndPoint &bdata,
        std::vector<std::pair<int, double>> &neibs,
        std::vector<double> &alpha,
        double &alpha0, double &beta,
        std::vector<Point> &bad
)
{
    double nux = bdata.n.x;
    double nuy = bdata.n.y;

    double omega = bdata.bc.a;
    double sigma = bdata.bc.b;
    double rho = bdata.bc.g;

    double e0 = 0;
    double e1 = 0;
    double e2 = 0;
    double e11 = 1;
    double e12 = 0;
    double e22 = 1;

    double f = 1;

    if (omega < 0) {
        omega = -omega;
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
        const Point &p = idata[j].first.p;
        double x = p.x - pc.x;
        double y = p.y - pc.y;
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
        logfile << "Failed to solve problem at boundary point " << pc << std::endl;
        for (const auto &jv : neibs) {
            int j = jv.first;
            const Point &p = idata[j].first.p;
            double x = p.x - pc.x;
            double y = p.y - pc.y;
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
        const Point &p = idata[neibs[j].first].first.p;
        double x = p.x - pc.x;
        double y = p.y - pc.y;
        double gj = 0.5 / (e11 + e22) * (x * x + y * y);

        alpha0 -= alp[j];
        alpha0 += alp[j] * gj * e0;
        beta += alp[j] * gj * f;
    }

    return true;
}

std::vector<Point> build_coeff(
    Interior inter, Boundary bnd, const Region *reg,
    std::vector<std::pair<IntPoint, Bucket<IntPoint> *> > &idata,
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
    std::vector<Point> bad;
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
        for (size_t j = i + 1; j < idata.size(); j++)
            if (idata[i].second->inall(idata[j].first.p)) { // same bucket
                double dd = idata[i].second->distance(idata[i].first.p, idata[j].first.p);

                const Point pa(idata[i].first.p);
                const Point pb(idata[j].first.p);

                for (double s = 0.05; s < 1; s += 0.1) {
                    Point mid(
                            s * pa.x + (1 - s) * pb.x,
                            s * pa.y + (1 - s) * pb.y
                        );
                    if (!reg->inside(mid))
                        dd = diam;
                }

                neibs[i].push_back(std::make_pair(j, dd));
                neibs[j].push_back(std::make_pair(i, dd));
            }
        bool isbnd = idata[i].first.bnd;
        if (!isbnd)
            continue;
        for (size_t j = 0; j < bdata.size(); j++)
            if (idata[i].second->incore(bdata[j].first.p)) {
                if (idata[i].second->distance(idata[i].first.p, bdata[j].first.p) < 1e-10)
                    bndidx[j] = i;
            }
    }
    for (size_t i = 0; i < idata.size(); i++) {
        std::sort(neibs[i].begin(), neibs[i].end(),
            [](const std::pair<int, double> &a, const std::pair<int, double> &b) -> bool {
                return a.second < b.second;
            }
        );
        if (neibs[i].size() < nneib)
            std::cerr << "Point " << idata[i].first.p << " has only " << neibs[i].size() << " neibs!" << std::endl;
        else
            neibs[i].erase(neibs[i].begin() + nneib, neibs[i].end());
    }

    int badint = 0;
    int badbnd = 0;

    for (size_t i = 0; i < neibs.size(); i++) {
        const Point &pc = idata[i].first.p;

//        std::cout << "p = " << pc<< " i = " << i << " Ns = " << neibs[i].size() << std::endl;
        if (!idata[i].first.bnd) {
            bool f = computeInternalEquation(pc, idata, neibs[i], alpha[i], alpha0[i], beta[i], bad);
            if (!f)
                badint++;
        }
    }

    for (size_t ib = 0; ib < bdata.size(); ib++) {
        size_t i = bndidx[ib];
        const Point &pc = bdata[ib].first.p;

//        std::cout << "p = " << pc << " i = " << i << " Ns = " << neibs[i].size() << std::endl;
        bool f = computeBoundaryEquation(pc, idata, bdata[ib].first, neibs[i], alpha[i], alpha0[i], beta[i], bad);
        if (!f)
            badbnd++;
    }

    std::cout << "Found " << badint << " bad internal points and " << badbnd << " bad boundary points" << std::endl;

    return bad;
}

int main(int argc, char **argv) {
	if (argc < 2) {
		std::cerr << "USAGE: driver <config file>" << std::endl;
		return 1;
	}

    logfile.open("bad_points.log");

	std::fstream f(argv[1], std::ios::in);
	Lexer lexer(f, argv[1]);

	const Region *config = 0;
	Parser parser(lexer, config);
//    parser.set_debug_level(10);

	if (parser.parse())
		return 1;

    std::cout << config->print() << std::endl;

    double hbord = 0.2;
    int lays = 3;
    double h = 0.2;

    const auto &bnds = config->bounds();
    Point ll = bnds.first;
    Point ur = bnds.second;

    ll.x -= h;
    ll.y -= h;
    ur.x -= h;
    ur.y -= h;

    Interior inter(ll, ur, 5 * h, 0.005);
    Boundary bnd(ll, ur, 5 * h, 0.005);

    mesh(config, bnd, inter, hbord, h, lays);

    std::vector<std::pair<IntPoint, Bucket<IntPoint> *>> idata;
    std::vector<int> bndidx;
    std::vector<std::vector<std::pair<int, double> > > neibs;
    std::vector<std::vector<double> > alpha;
    std::vector<double> alpha0;
    std::vector<double> beta;

    int maxiter = 10;

    for (int iter = 0; ; iter++) {
        const auto bad = build_coeff(inter, bnd, config,
                idata, bndidx, neibs, alpha, alpha0, beta
            );
        std::cout << "From " << idata.size() << " points "<< bad.size() << " are marked as bad" << std::endl;

        std::ofstream f("neibs." + std::to_string(iter));
        for (size_t i = 0; i < neibs.size(); i++)
            if (!idata[i].first.bnd)
                for (const auto &jv : neibs[i]) {
                    int j = jv.first;
                    f << idata[i].first.p.x << " " << idata[i].first.p.y << std::endl;
                    f << idata[j].first.p.x << " " << idata[j].first.p.y << std::endl << std::endl;
                }
        f.close();

        std::ofstream ff("bnds." + std::to_string(iter));
        for (size_t i = 0; i < neibs.size(); i++)
            if (idata[i].first.bnd)
                for (const auto &jv : neibs[i]) {
                    int j = jv.first;
                    ff << idata[i].first.p.x << " " << idata[i].first.p.y << std::endl;
                    ff << idata[j].first.p.x << " " << idata[j].first.p.y << std::endl << std::endl;
                }
        ff.close();

        std::ofstream g("bad." + std::to_string(iter));
        for (size_t i = 0; i < bad.size(); i++)
            g << bad[i].x << " " << bad[i].y << std::endl;
        g.close();

        if (bad.empty())
            break;

        for (const auto &p : bad) {
            inter.remove(p);
            bnd.remove(p);
        }

        if (iter > maxiter) {
            std::cerr << "Could not remove all bad points in " << maxiter << " iterations. Try different stepsize" << std::endl;
            abort();
        }
    }

	return 0;
}
