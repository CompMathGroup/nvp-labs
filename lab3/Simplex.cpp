#include "Simplex.h"

#include <cassert>
#include <Eigen/Core>
#include <Eigen/LU>
#include <iostream>
#include <limits>

#include <stdio.h>

using namespace std;
using namespace Eigen;

namespace simplex {

static const double tol = 1e-12;

template <typename DA, typename Db, typename Dc>
void print_program(
        const MatrixBase<DA> &A,
        const MatrixBase<Db> &b,
        const MatrixBase<Dc> &c,
        const vector<int> &basis,
        const double p
        )
{
    int m = b.size();
    int n = c.size();
    printf("bvar     obj   |");
    for (int i = 0; i < n / 2; i++)
        printf("          ");
    printf("A");
    for (int i = n / 2; i < n - 1; i++)
        printf("          ");
    printf("         |      b\n");
    for (int i = 0; i < m; i++) {
        printf("%2d |           |", basis[i]);
        for (int j = 0; j < n; j++)
            printf("% 8.2e ", A(i, j));
        printf("| % 8.2e\n", b[i]);
    }
    printf("---+-----------+");
    for (int i = 0; i < n; i++)
        printf("----------");
    printf("+-----------\n");
    printf(" \xce\x94 | % 8.2e |", p);
    for (int j = 0; j < n; j++)
        printf("% 8.2e ", c[j]);
    printf("\n");
    printf("\n");
}

/* Removes linearly dependent constraints
 * Negates some equations so every b is positive
 * Return true if can proceed (! (es & Inconsistent))
 * */
bool simplifyConstraints(MatrixXd &A, VectorXd &b, ExtraStatus &es) {
    int m = A.rows();
    FullPivLU<MatrixXd> lu(A);

    int rank = lu.rank();
    es = ExtraStatus::Normal;
    if (rank < m)
        es = es | ExtraStatus::IncompleteRank;

    const auto &LU = lu.matrixLU();
    const auto &L = LU.block(0, 0, m, m).triangularView<UnitLower>();
    const auto &P = lu.permutationP();
    const auto &Q = lu.permutationQ();

    A = LU.triangularView<Upper>();
    A = A * Q.transpose();
    b = L.solve(P * b);
    double bthresh = 10 * lu.threshold() * b.norm();

    for (int i = rank; i < m; i++)
        if (std::abs(b[i]) > bthresh) {
            es = es | ExtraStatus::Inconsistent;
            return false;
        }

    A.conservativeResize(rank, NoChange);
    b.conservativeResize(rank, NoChange);
    m = rank;

    for (int i = 0; i < m; i++)
        if (b[i] < -bthresh) {
            A.row(i) *= -1;
            b[i] = -b[i];
        }

    return true;
}

Status simplexMethod(MatrixXd &A, VectorXd &b, RowVectorXd &c, vector<int> &basis, double &obj, bool debug) {
    const double cthresh = tol * c.norm();
    const double athresh = tol * A.norm();

    int m = A.rows();
    int n = A.cols();

    for (int i = 0; i < m; i++) {
        int j = basis[i];
        c -= A.row(i) * c[j];
    }

    const int maxiters = 1000;
    for (int iter = 0; iter < maxiters; iter++) {
        if (debug)
            print_program(A, b, c, basis, obj);
        int ev = 0;
        for (ev = 0; ev < n; ev++)
            if (c[ev] < -cthresh)
                break;

        if (ev == n) /* Cannot optimize further -> solution */
            return Status::Solved;

        int lvr = -1;
        double minratio;
        for (int i = 0; i < m; i++)
            if (A(i, ev) > athresh) {
                double ratio = b[i] / A(i, ev);
                if (lvr < 0 || ratio < minratio) {
                    lvr = i;
                    minratio = ratio;
                }
            }

        if (lvr == -1)
            return Status::Unbounded;

        basis[lvr] = ev;
        double p = A(lvr, ev);

        A.row(lvr) /= p;
        b[lvr] /= p;

        for (int i = 0; i < m; i++) {
            if (i == lvr)
                continue;

            double p2 = A(i, ev);
            A.row(i) -= p2 * A.row(lvr);
            b[i] -= p2 * b[lvr];
        }

        double p2 = c[ev];
        c -= p2 * A.row(lvr);
        obj += p2 * b[lvr];
    }

    return Status::FailedToStop;
}

/* On exit A(:, basis) = E
 *
 * */
bool findFeasiblePoint(MatrixXd &A, VectorXd &b, vector<int> &basis, bool debug) {
    int m = A.rows();
    int n = A.cols();
    MatrixXd AE(m, m + n);
    AE.block(0, 0, m, n) = A;
    AE.block(0, n, m, m) = MatrixXd::Identity(m, m);
    vector<int> extbasis(m);
    for (int i = 0; i < m; i++)
        extbasis[i] = i + n;
    RowVectorXd cc(m + n);
    for (int i = 0; i < n; i++)
        cc[i] = 0;
    for (int i = n; i < n + m; i++)
        cc[i] = 1;
    double obj = cc.block(0, n, 1, m).transpose().dot(b);

    /* The problem is consistent and feasible. Also it is bounded. */
    Status phase1 = simplexMethod(AE, b, cc, extbasis, obj, debug);

    if (phase1 != Status::Solved) {
        cerr << "Phase 1 failed!" << endl;
        return false;
    }

    if (obj > 10 * m * tol)
        return false;

//    bool spurious = false;
    A = AE.block(0, 0, m, n);

    for (int i = 0; i < m; i++) {
        if (extbasis[i] >= n) {
            int jmax = -1;
            double maxv = 0;
            for (int j = 0; j < n; j++) {
                if (std::abs(A(i, j)) > maxv) {
                    jmax = j;
                    maxv = std::abs(A(i, j));
                }
            }
            double piv = A(i, jmax);
            A.row(i) /= piv;
            b[i] /= piv;
            for (int ii = 0; ii < m; ii++)
                if (i != ii) {
                    double pp = A(ii, jmax);
                    A.row(ii) -= pp * A.row(i);
                    b[ii] -= pp * b[i];
                }
            extbasis[i] = jmax;
//            spurious = true;
        }
        basis[i] = extbasis[i];
    }
/*
    if (spurious)
        cout << "Spurious variables at optimum" << endl;
*/
    return true;
}

Status solveProblem(
        const vector<double> &_A,
        const vector<double> &_b,
        const vector<double> &_c,
        vector<int> &basis,
        vector<double> &_x,
        double *p,
        ExtraStatus *pes,
        bool debug
    )
{
    size_t m = _b.size();
    size_t n = _c.size();
    assert(m * n == _A.size());

    MatrixXd A(m, n);
    VectorXd b(m);
    RowVectorXd c(n);

    copy(_A.begin(), _A.end(), A.data());
    copy(_b.begin(), _b.end(), b.data());
    copy(_c.begin(), _c.end(), c.data());

    ExtraStatus es;
    simplifyConstraints(A, b, es);
    m = A.rows();

    if (pes)
        *pes = es;

    if (es & ExtraStatus::Inconsistent) {
        return Status::Infeasible;
    }

    basis.resize(m);
    bool ret = findFeasiblePoint(A, b, basis, debug);

    if (!ret)
        return Status::Infeasible;

    double obj = 0;
    for (size_t i = 0; i < m; i++)
        obj += c[basis[i]] * b[i];
    Status phase2 = simplexMethod(A, b, c, basis, obj, debug);

    _x.resize(m);
    copy(b.data(), b.data() + m, _x.begin());

    if (p)
        *p = obj;

    return phase2;
}

}
/*
int main() {
    vector<double> A = {
        -1, 2, 1, 2,
        0, 3, 3, 6,
        4, 3, 1, 8,
        3, 5, 2, 10,
        0, -1, 1, 0
    };
    vector<double> b = {2, 3, 2, 7};
    vector<double> c = {-1, -3, -2, -4, -1};
    vector<int> basis;
    vector<double> x;

    simplex::ExtraStatus es;
    double obj;

    simplex::Status s = simplex::solveProblem(A, b, c, basis, x, &obj, &es);

    cout << "Status: " << static_cast<int>(s) << ", ExtraStatus: " << static_cast<int>(s) << ", p = " << obj << endl;

    for (size_t i = 0; i < x.size(); i++)
        cout << "x[" << basis[i] << "] = " << x[i] << endl;
}
*/
