#include <Eigen/Core>
#include <Eigen/LU>
#include <cmath>
#include <iostream>
#include <fstream>

#include <fenv.h>

template<int n>
using vec = Eigen::Matrix<double, n, 1>;

template<int n>
using mat = Eigen::Matrix<double, n, n>;

inline double sqr(double x) {
    return x*x;
}

template<int n>
void eno5(
        const vec<n> U[5],
        const mat<n> Omega[4],
        const mat<n> invOmega[4],
        vec<n> &Uplus,
        vec<n> &Uminus
    )
{
    vec<n> dW[4];
    vec<n> AdWp[4];
    vec<n> AdWm[4];

    for (int i = 0; i < 4; i++) {
        dW[i] = Omega[i] * (U[i+1] - U[i]);
    }

    for (int i = 0; i < n; i++) {
        double dmm = dW[0][i];
        double dm  = dW[1][i];
        double dp  = dW[2][i];
        double dpp = dW[3][i];

        double b1 = (3 * sqr(3*dm - dmm) + 13 * sqr(dm - dmm)) / 48;
        double b2 = (3 * sqr(dm + dp) + 13 * sqr(dp - dm)) / 48;
        double b3 = (3 * sqr(3*dp - dpp) + 13 * sqr(dpp - dp)) / 48;

        int best = 0;
        double bmin = b1;
        if (b2 < bmin) {
            best = 1;
            bmin = b2;
        }
        if (b3 < bmin) {
            best = 2;
            bmin = b3;
        }

        double a[3] = {0, 0, 0};
        a[best] = 1;

        AdWp[0][i] = -a[0] / 3 * dmm;
        AdWp[1][i] = (5*a[0]+a[1]) / 6 * dm;
        AdWp[2][i] = (a[1]+2*a[2]) / 3 * dp;
        AdWp[3][i] = -a[2] / 6 * dpp;

        AdWm[0][i] = a[0] / 6 * dmm;
        AdWm[1][i] = -(2*a[0]+a[1]) / 3 * dm;
        AdWm[2][i] = -(a[1]+5*a[2]) / 6 * dp;
        AdWm[3][i] = a[2] / 3 * dpp;
    }

    Uplus  = U[2];
    Uminus = U[2];

    for (int i = 0; i < 4; i++) {
        Uplus  += invOmega[i] * AdWp[i];
        Uminus += invOmega[i] * AdWm[i];
    }
}

const double GAMMA = 1.4;

vec<3> from_rup(double rho, double u, double p) {
    return vec<3>(rho, rho * u, 0.5 * rho * u * u + p / (GAMMA - 1));
}

mat<3> comp_Omega(const vec<3> &U, vec<3> &lam) {
    const double u = U[1] / U[0];
    const double K = .5 * U[1] * u;
    double p = (GAMMA - 1) * (U[2] - K);
    double c = std::sqrt(GAMMA * p / U[0]);

    double u2 = u * u;
    double c2 = c * c;
    double zm = (GAMMA - 1) / (c2 + (GAMMA - 1) * u * (0.5 * u - c));
    double zp = (GAMMA - 1) / (c2 + (GAMMA - 1) * u * (0.5 * u + c));

    mat<3> W;

    W(0, 0) = 1;        W(0, 1) = zm;           W(0, 2) = zp;
    W(1, 0) = u;        W(1, 1) = zm * (u - c); W(1, 2) = zp * (u + c);
    W(2, 0) = .5 * u2;  W(2, 1) = 1;            W(2, 2) = 1;

    lam[0] = u;
    lam[1] = u - c;
    lam[2] = u + c;

    return W;
}

double comp_amax(const vec<3> &U) {
    const double u = U[1] / U[0];
    const double K = .5 * U[1] * u;
    const double p = (GAMMA - 1) * (U[2] - K);
    const double c = std::sqrt(GAMMA * p / U[0]);

    return c + std::abs(u);
}

vec<3> flux(const vec<3> &U) {
    const double u = U[1] / U[0];
    const double K = .5 * U[1] * u;
    const double p = (GAMMA - 1) * (U[2] - K);

    return vec<3>(U[1], p + U[1] * u, (U[2] + p) * u);
}

vec<3> riemann(const vec<3> &UL, const vec<3> &UR, const mat<3> &W, const vec<3> &lam, const mat<3> &iW) {
    auto dW = W * (UL - UR);
    const vec<3> alam(std::abs(lam[0]), std::abs(lam[1]), std::abs(lam[2]));

    return 0.5 * (flux(UL) + flux(UR)) + 0.5 * iW * (alam.cwiseProduct(dW));
}

void euler_step(const double dt_h, const int M, const vec<3> U[], vec<3> Unext[]) {
    vec<3> Um[M], Up[M];
    vec<3> lam[M-1];
    mat<3> Omega[M-1], invOmega[M-1];
    vec<3> F[M-1];

    for (int j = 0; j < M-1; j++) {
        vec<3> Uavg = 0.5 * (U[j] + U[j+1]);

        invOmega[j] = comp_Omega(Uavg, lam[j]);
        Omega[j] = invOmega[j].inverse();
    }

    for (int j = 0; j < 2; j++)
        Um[j] = Up[j] = U[j];
    for (int j = M-2; j < M; j++)
        Um[j] = Up[j] = U[j];
    for (int j = 2; j < M-2; j++)
        eno5(&U[j-2], &Omega[j-2], &invOmega[j-2], Up[j], Um[j]);
//        Um[j] = Up[j] = U[j];

    for (int j = 0; j < M-1; j++)
        F[j] = riemann(Up[j], Um[j+1], Omega[j], lam[j], invOmega[j]);

    Unext[0]   = U[0];
    Unext[M-1] = U[M-1];

    for (int j = 1; j < M-1; j++)
        Unext[j] = U[j] + dt_h * (F[j-1] - F[j]);
}

int main() {
    feenableexcept(FE_INVALID | FE_OVERFLOW | FE_DIVBYZERO);

    const int M = 500;
    const double h = 1. / M;
    const double C = 0.4;

    const int N = 400;

    vec<3> U[M];
    vec<3> U1[M];
    vec<3> U2[M];

    for (int i = 0; i < M; i++) {
        double xi = (i + .5) * h - .5;
        if (xi < 0)
            U[i] = from_rup(1, 0.75, 1);
        else
            U[i] = from_rup(0.125, 0, 0.1);
    }

    for (int n = 0; n < N; n++) {
        double amax = 0;
        for (int j = 0; j < M; j++)
            amax = std::max(amax, comp_amax(U[j]));
        const double dt = C * h / amax;

        euler_step(dt / h, M, U, U1);
/*
        for (int j = 0; j < M; j++)
            U[j] = U1[j];
*/
        euler_step(dt / h, M, U1, U2);

        for (int j = 0; j < M; j++)
            U[j] = 0.5 * (U1[j] + U2[j]);

    }

    std::ofstream f("res.csv");
    f << "x,rho\n";
    for (int i = 0; i < M; i++) {
        double xi = (i + .5) * h - .5;
        f << xi << "," << U[i][0] << "\n";
    }

    return 0;
}
