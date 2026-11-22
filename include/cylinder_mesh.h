#pragma once
// Polar structured mesh for the annulus a < r < R.
// Nodes: (Nr+1) radial layers x Nt angular divisions, periodic in theta.
// Each quad is split into 2 triangles.
#include <algorithm>
#include <cmath>
#include <map>
#include <vector>

struct Node {
    double r = 0.0, theta = 0.0, x = 0.0, y = 0.0;
};

struct Element {
    int n[3] = {0, 0, 0};  // node ids, CCW
    double area = 0.0;
    double grad[3][2] = {};  // grad of P1 basis, [node][x/y]
};

struct Edge {
    int n0 = -1, n1 = -1;
    int elem0 = -1, elem1 = -1;  // second is -1 on the boundary
    double length = 0.0;
    double mx = 0.0, my = 0.0;  // midpoint
    bool is_inner = false;      // r == a, PEC
    bool is_outer = false;      // r == R, Robin
};

struct Mesh {
    double a = 1.0, R = 3.0;
    int Nr = 4, Nt = 16;
    std::vector<Node> nodes;
    std::vector<Element> elems;
    std::vector<Edge> edges;
    // per element, its 3 edge ids, in the same order as its nodes:
    // edge i is opposite node i, i.e. between the other two nodes.
    std::vector<int> elem_edges;  // size = elems.size() * 3
};

inline Mesh build_cylinder_mesh(double a, double R, int Nr, int Nt) {
    Mesh m;
    m.a = a; m.R = R; m.Nr = Nr; m.Nt = Nt;
    const double pi = std::acos(-1.0);

    auto node_id = [&](int ir, int it) {
        return ir * Nt + ((it % Nt + Nt) % Nt);
    };

    for (int ir = 0; ir <= Nr; ++ir) {
        double r = a + (R - a) * ir / Nr;
        for (int it = 0; it < Nt; ++it) {
            double th = 2.0 * pi * it / Nt;
            Node nd;
            nd.r = r; nd.theta = th;
            nd.x = r * std::cos(th); nd.y = r * std::sin(th);
            m.nodes.push_back(nd);
        }
    }

    auto add_tri = [&](int n0, int n1, int n2) {
        Element e;
        e.n[0] = n0; e.n[1] = n1; e.n[2] = n2;
        const Node &p0 = m.nodes[n0], &p1 = m.nodes[n1], &p2 = m.nodes[n2];
        double x0 = p0.x, y0 = p0.y, x1 = p1.x, y1 = p1.y, x2 = p2.x, y2 = p2.y;
        double det = (x1 - x0) * (y2 - y0) - (x2 - x0) * (y1 - y0);
        e.area = 0.5 * std::abs(det);
        // Gradients of P1 basis functions.
        double inv2A = 1.0 / det;  // signed, keeps orientation consistent
        e.grad[0][0] = (y1 - y2) * inv2A;
        e.grad[0][1] = (x2 - x1) * inv2A;
        e.grad[1][0] = (y2 - y0) * inv2A;
        e.grad[1][1] = (x0 - x2) * inv2A;
        e.grad[2][0] = (y0 - y1) * inv2A;
        e.grad[2][1] = (x1 - x0) * inv2A;
        m.elems.push_back(e);
    };

    for (int ir = 0; ir < Nr; ++ir) {
        for (int it = 0; it < Nt; ++it) {
            int it1 = (it + 1) % Nt;
            int n00 = node_id(ir, it), n01 = node_id(ir, it1);
            int n10 = node_id(ir + 1, it), n11 = node_id(ir + 1, it1);
            add_tri(n00, n10, n11);
            add_tri(n00, n11, n01);
        }
    }

    // Build edges from elements.
    std::map<std::pair<int, int>, int> edge_map;
    for (int ei = 0; ei < (int)m.elems.size(); ++ei) {
        const Element &e = m.elems[ei];
        for (int k = 0; k < 3; ++k) {
            int a0 = e.n[(k + 1) % 3], a1 = e.n[(k + 2) % 3];
            auto key = std::minmax(a0, a1);
            auto itf = edge_map.find(key);
            int eid;
            if (itf == edge_map.end()) {
                Edge ed;
                ed.n0 = key.first; ed.n1 = key.second;
                ed.elem0 = ei;
                const Node &p0 = m.nodes[ed.n0], &p1 = m.nodes[ed.n1];
                ed.length = std::hypot(p1.x - p0.x, p1.y - p0.y);
                ed.mx = 0.5 * (p0.x + p1.x);
                ed.my = 0.5 * (p0.y + p1.y);
                // Midpoint radius is inside the chord, so test the nodes, not the midpoint.
                ed.is_inner = std::abs(p0.r - a) < 1e-9 && std::abs(p1.r - a) < 1e-9;
                ed.is_outer = std::abs(p0.r - R) < 1e-9 && std::abs(p1.r - R) < 1e-9;
                eid = (int)m.edges.size();
                m.edges.push_back(ed);
                edge_map[key] = eid;
            } else {
                eid = itf->second;
                m.edges[eid].elem1 = ei;
            }
            // elem_edges is filled after the loop in element order below.
            (void)eid;
        }
    }
    m.elem_edges.assign(m.elems.size() * 3, -1);
    for (int ei = 0; ei < (int)m.elems.size(); ++ei) {
        const Element &e = m.elems[ei];
        for (int k = 0; k < 3; ++k) {
            int a0 = e.n[(k + 1) % 3], a1 = e.n[(k + 2) % 3];
            auto key = std::minmax(a0, a1);
            m.elem_edges[ei * 3 + k] = edge_map[key];
        }
    }
    return m;
}
