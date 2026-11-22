#pragma once
// Simplified HDG for the TM Helmholtz problem.
// Per element (P1 Ez, P0 trace lambda on each edge):
//   A = K - k^2 M + tau * E
//   B[i][e] = tau * int_edge phi_i ds
//   C[e][e] = tau * |edge|
// Static condensation: S = C - B^T A^{-1} B, solve for lambda, recover Ez.
// This is a readable demo formulation, not a validated solver.
// TODO: double-trace HDG (tangential E and H traces as hybrid unknowns).
#include <complex>
#include <vector>

using cx = std::complex<double>;

// Invert a 3x3 complex matrix by adjugate. Returns false if singular.
inline bool invert3(const cx a[3][3], cx inv[3][3]) {
    cx det = a[0][0] * (a[1][1] * a[2][2] - a[1][2] * a[2][1])
           - a[0][1] * (a[1][0] * a[2][2] - a[1][2] * a[2][0])
           + a[0][2] * (a[1][0] * a[2][1] - a[1][1] * a[2][0]);
    if (std::abs(det) < 1e-14) return false;
    inv[0][0] = (a[1][1] * a[2][2] - a[1][2] * a[2][1]) / det;
    inv[0][1] = (a[0][2] * a[2][1] - a[0][1] * a[2][2]) / det;
    inv[0][2] = (a[0][1] * a[1][2] - a[0][2] * a[1][1]) / det;
    inv[1][0] = (a[1][2] * a[2][0] - a[1][0] * a[2][2]) / det;
    inv[1][1] = (a[0][0] * a[2][2] - a[0][2] * a[2][0]) / det;
    inv[1][2] = (a[0][2] * a[1][0] - a[0][0] * a[1][2]) / det;
    inv[2][0] = (a[1][0] * a[2][1] - a[1][1] * a[2][0]) / det;
    inv[2][1] = (a[0][1] * a[2][0] - a[0][0] * a[2][1]) / det;
    inv[2][2] = (a[0][0] * a[1][1] - a[0][1] * a[1][0]) / det;
    return true;
}

// Dense complex solve Ax = b by Gaussian elimination with partial pivoting.
// A is n x n row-major, overwritten. Returns false if singular.
inline bool solve_dense(std::vector<cx> &A, std::vector<cx> &b, int n) {
    for (int col = 0; col < n; ++col) {
        int piv = col;
        double best = std::abs(A[col * n + col]);
        for (int r = col + 1; r < n; ++r) {
            double v = std::abs(A[r * n + col]);
            if (v > best) { best = v; piv = r; }
        }
        if (best < 1e-14) return false;
        if (piv != col) {
            for (int c = col; c < n; ++c) std::swap(A[piv * n + c], A[col * n + c]);
            std::swap(b[piv], b[col]);
        }
        cx diag = A[col * n + col];
        for (int r = col + 1; r < n; ++r) {
            cx f = A[r * n + col] / diag;
            if (f == cx(0.0, 0.0)) continue;
            for (int c = col; c < n; ++c) A[r * n + c] -= f * A[col * n + c];
            b[r] -= f * b[col];
        }
    }
    for (int r = n - 1; r >= 0; --r) {
        cx s = b[r];
        for (int c = r + 1; c < n; ++c) s -= A[r * n + c] * b[c];
        b[r] = s / A[r * n + r];
    }
    return true;
}
