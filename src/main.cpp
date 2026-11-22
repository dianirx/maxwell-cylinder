// hdg-maxwell-cylinder: minimal HDG demo for TM scattering by a cylinder.
// Simplified P1/P0 HDG, static condensation, dense trace solve.
// Demo only - not a validated solver, no results claimed.
#include <cmath>
#include <complex>
#include <fstream>
#include <iostream>
#include <vector>

#include "cylinder_mesh.h"
#include "hdg_tm.h"

using cx = std::complex<double>;

int main() {
    const double a = 1.0;   // inner cylinder radius
    const double R = 3.0;   // outer boundary radius
    const int Nr = 4, Nt = 16;
    const double freq = 1.0e9;         // 1 GHz
    const double c0 = 3.0e8;
    const double pi = std::acos(-1.0);
    const double k = 2.0 * pi * freq / c0;
    const cx tau(k, 0.0);
    const cx I(0.0, 1.0);

    Mesh mesh = build_cylinder_mesh(a, R, Nr, Nt);
    std::cout << "HDG TM cylinder demo\n";
    std::cout << "  a=" << a << " R=" << R << " k=" << k << " Nr=" << Nr << " Nt=" << Nt << "\n";
    std::cout << "  nodes=" << mesh.nodes.size()
              << " elems=" << mesh.elems.size()
              << " edges=" << mesh.edges.size() << "\n";

    // Unknowns: lambda on every edge except inner PEC edges (lambda = 0 there).
    std::vector<int> edge_unknown(mesh.edges.size(), -1);
    int n_unknown = 0;
    for (int e = 0; e < (int)mesh.edges.size(); ++e) {
        if (!mesh.edges[e].is_inner) edge_unknown[e] = n_unknown++;
    }
    std::cout << "  trace unknowns=" << n_unknown << "\n";

    std::vector<cx> S((size_t)n_unknown * n_unknown, cx(0.0, 0.0));
    std::vector<cx> rhs(n_unknown, cx(0.0, 0.0));
    std::vector<cx> elem_Ez(mesh.elems.size() * 3, cx(0.0, 0.0));

    // Assembly: one element at a time, Schur complement.
    for (int ei = 0; ei < (int)mesh.elems.size(); ++ei) {
        const Element &el = mesh.elems[ei];
        double area = el.area;

        // P1 mass: area/12 * [[2,1,1],[1,2,1],[1,1,2]]
        // P1 stiffness: area * grad_i . grad_j
        cx A[3][3], B[3][3], C[3][3];
        for (int i = 0; i < 3; ++i) {
            for (int j = 0; j < 3; ++j) {
                double stiff = area * (el.grad[i][0] * el.grad[j][0] +
                                       el.grad[i][1] * el.grad[j][1]);
                double mass = area / 12.0 * ((i == j) ? 2.0 : 1.0);
                A[i][j] = stiff - k * k * mass;
                B[i][j] = cx(0.0, 0.0);
                C[i][j] = cx(0.0, 0.0);
            }
        }
        // Edge contributions. Edge j of the element is opposite node j.
        // int_edge phi_i ds = |edge|/2 if node i is an endpoint, else 0.
        // int_edge phi_i phi_j ds = |edge|/3 (i==j, endpoint), |edge|/6 (both endpoints).
        cx g[3] = {cx(0.0, 0.0), cx(0.0, 0.0), cx(0.0, 0.0)};
        for (int ej = 0; ej < 3; ++ej) {
            int eid = mesh.elem_edges[ei * 3 + ej];
            const Edge &ed = mesh.edges[eid];
            double L = ed.length;
            int nA = el.n[(ej + 1) % 3], nB = el.n[(ej + 2) % 3];
            (void)nA; (void)nB;
            for (int i = 0; i < 3; ++i) {
                bool endpoint = (i == (ej + 1) % 3) || (i == (ej + 2) % 3);
                if (endpoint) B[i][ej] += tau * (L / 2.0);
                for (int j = 0; j < 3; ++j) {
                    bool j_endpoint = (j == (ej + 1) % 3) || (j == (ej + 2) % 3);
                    if (endpoint && j_endpoint) {
                        double v = (i == j) ? L / 3.0 : L / 6.0;
                        A[i][j] += tau * v;
                    }
                }
            }
            C[ej][ej] += tau * L;
            if (ed.is_outer) {
                // Robin: + i k |edge| on C, rhs 2 i k |edge| Ez_inc(midpoint).
                cx ez_inc = std::exp(I * k * ed.mx);  // exp(i k x), x = mx
                C[ej][ej] += I * k * L;
                g[ej] += cx(2.0, 0.0) * I * k * L * ez_inc;
            }
        }

        cx Ainv[3][3];
        if (!invert3(A, Ainv)) {
            std::cerr << "singular local matrix on element " << ei << "\n";
            return 1;
        }
        // S_e = C - B^T A^{-1} B, rhs_e = g (f = 0)
        cx Se[3][3];
        for (int i = 0; i < 3; ++i) {
            for (int j = 0; j < 3; ++j) {
                cx s = C[i][j];
                for (int p = 0; p < 3; ++p)
                    for (int q = 0; q < 3; ++q)
                        s -= B[p][i] * Ainv[p][q] * B[q][j];
                Se[i][j] = s;
            }
        }
        int eids[3] = {mesh.elem_edges[ei * 3 + 0], mesh.elem_edges[ei * 3 + 1],
                       mesh.elem_edges[ei * 3 + 2]};
        for (int i = 0; i < 3; ++i) {
            int ui = edge_unknown[eids[i]];
            if (ui < 0) continue;  // PEC edge, lambda = 0, no unknown
            rhs[ui] += g[i];
            for (int j = 0; j < 3; ++j) {
                int uj = edge_unknown[eids[j]];
                if (uj < 0) continue;
                S[(size_t)ui * n_unknown + uj] += Se[i][j];
            }
        }
        // Store Ainv and B for the recovery pass (recompute cheaply below instead).
        (void)Ainv;
    }

    if (!solve_dense(S, rhs, n_unknown)) {
        std::cerr << "global trace solve failed (singular)\n";
        return 1;
    }
    std::vector<cx> lambda(mesh.edges.size(), cx(0.0, 0.0));
    for (int e = 0; e < (int)mesh.edges.size(); ++e)
        if (edge_unknown[e] >= 0) lambda[e] = rhs[edge_unknown[e]];

    // Recovery pass: Ez per element, u = A^{-1} (f - B lambda), f = 0.
    double max_mag = 0.0;
    std::ofstream out("field.csv");
    out << "x,y,Ez_real,Ez_imag,abs_Ez\n";
    for (int ei = 0; ei < (int)mesh.elems.size(); ++ei) {
        const Element &el = mesh.elems[ei];
        double area = el.area;
        cx A[3][3], B[3][3];
        for (int i = 0; i < 3; ++i) {
            for (int j = 0; j < 3; ++j) {
                double stiff = area * (el.grad[i][0] * el.grad[j][0] +
                                       el.grad[i][1] * el.grad[j][1]);
                double mass = area / 12.0 * ((i == j) ? 2.0 : 1.0);
                A[i][j] = stiff - k * k * mass;
                B[i][j] = cx(0.0, 0.0);
            }
        }
        for (int ej = 0; ej < 3; ++ej) {
            int eid = mesh.elem_edges[ei * 3 + ej];
            const Edge &ed = mesh.edges[eid];
            double L = ed.length;
            for (int i = 0; i < 3; ++i) {
                bool endpoint = (i == (ej + 1) % 3) || (i == (ej + 2) % 3);
                if (endpoint) B[i][ej] += tau * (L / 2.0);
                for (int j = 0; j < 3; ++j) {
                    bool j_endpoint = (j == (ej + 1) % 3) || (j == (ej + 2) % 3);
                    if (endpoint && j_endpoint) {
                        double v = (i == j) ? L / 3.0 : L / 6.0;
                        A[i][j] += tau * v;
                    }
                }
            }
        }
        cx Ainv[3][3];
        if (!invert3(A, Ainv)) return 1;
        cx lam[3] = {lambda[mesh.elem_edges[ei * 3 + 0]],
                     lambda[mesh.elem_edges[ei * 3 + 1]],
                     lambda[mesh.elem_edges[ei * 3 + 2]]};
        cx cx_cent(0.0, 0.0), cy_cent(0.0, 0.0);
        for (int i = 0; i < 3; ++i) {
            cx u = cx(0.0, 0.0);
            for (int p = 0; p < 3; ++p)
                for (int q = 0; q < 3; ++q)
                    u -= Ainv[i][p] * B[p][q] * lam[q];
            elem_Ez[ei * 3 + i] = u;
            const Node &nd = mesh.nodes[el.n[i]];
            cx_cent += nd.x / 3.0;
            cy_cent += nd.y / 3.0;
        }
        cx ez_avg = (elem_Ez[ei * 3 + 0] + elem_Ez[ei * 3 + 1] + elem_Ez[ei * 3 + 2]) / 3.0;
        max_mag = std::max(max_mag, std::abs(ez_avg));
        out << std::real(cx_cent) << "," << std::real(cy_cent) << ","
            << std::real(ez_avg) << "," << std::imag(ez_avg) << ","
            << std::abs(ez_avg) << "\n";
    }
    out.close();

    // Residual check of the trace system (re-assemble is skipped; report magnitude).
    std::cout << "  solved, max |Ez| (centroid avg) = " << max_mag << "\n";
    std::cout << "  sample Ez at element 0: " << elem_Ez[0] << ", " << elem_Ez[1]
              << ", " << elem_Ez[2] << "\n";
    std::cout << "  wrote field.csv (" << mesh.elems.size() << " rows)\n";
    std::cout << "Done. Demo only, results not validated.\n";
    return 0;
}
