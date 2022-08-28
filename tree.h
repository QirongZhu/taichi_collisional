#include "common.h"
#include "kernel.h"

#include <CGAL/Simple_cartesian.h>
#include <CGAL/Min_sphere_of_points_d_traits_3.h>
#include <CGAL/Min_sphere_of_spheres_d.h>

typedef CGAL::Simple_cartesian<double> Kp;
typedef CGAL::Min_sphere_of_points_d_traits_3<Kp, double> Traitsp;
typedef CGAL::Min_sphere_of_spheres_d<Traitsp> Min_spherep;
typedef Kp::Point_3 Pointp;

#include <CGAL/Min_sphere_of_spheres_d.h>
#include <CGAL/Cartesian.h>
#include <CGAL/Min_sphere_of_spheres_d.h>

typedef double FT;
typedef CGAL::Cartesian<FT> K;
typedef CGAL::Min_sphere_of_spheres_d_traits_3<K, FT> Traits;
typedef CGAL::Min_sphere_of_spheres_d<Traits> Min_sphere;
typedef K::Point_3 Point;
typedef Traits::Sphere Sphere;

namespace FMM
{

    class Tree
    {
    public:

        static void initKernel()
        {
            for (int n = 0; n <= EXPANSION; n++)
            {

                factorial_table[0] = 1.0;

                real_t inverse_factorial_table[4 * EXPANSION + 1];
                inverse_factorial_table[0] = 1.0;

                for (int i = 1; i <= EXPANSION * 2; i++)
                {
                    factorial_table[i] = (real_t)i * factorial_table[i - 1];
                    inverse_factorial_table[i] = 1 / factorial_table[i];
                }

                for (int i = 0; i <= EXPANSION; i++)
                {
                    for (int j = 0; j <= EXPANSION; j++)
                    {
                        factorial_coef[i][j] = factorial_table[i + j] * factorial_table[i - j];
                        factorial_coef_inv[i][j] = 1.0 / factorial_coef[i][j];
                        if (j <= i)
                            combinator_coef[i][j] =
                                factorial_table[i] * inverse_factorial_table[j] * inverse_factorial_table[i - j];
                    }
                }

                for (int n = 0; n <= EXPANSION; n++)
                {
                    int index_start = n * n + n;
                    factorial_coef_oned[index_start] = factorial_coef[n][0];
                    factorial_coef_inv_oned[index_start] = factorial_coef_inv[n][0];
                    for (int m = 1; m <= n; m++)
                    {
                        factorial_coef_oned[index_start + m] = factorial_coef[n][m];
                        factorial_coef_oned[index_start - m] = factorial_coef[n][m];
                        factorial_coef_inv_oned[index_start + m] = factorial_coef_inv[n][m];
                        factorial_coef_inv_oned[index_start - m] = factorial_coef_inv[n][m];
                    }
                }
            }
        }

        int findSplitAxis(Biter begin, Biter end, double &split_pos)
        {
            double mx = 0, my = 0, mz = 0, mtot = 0;
            double minx = 1e30, miny = 1e30, minz = 1e30;
            double maxx = -1e30, maxy = -1e30, maxz = -1e30;
            auto it = begin;

            for (it = begin; it != end; it++)
            {
                auto m = (*it).m;
                auto x = (*it).X[0];
                auto y = (*it).X[1];
                auto z = (*it).X[2];

                mtot += m;

                mx += x * m;
                my += y * m;
                mz += z * m;

                minx = x < minx ? x : minx;
                maxx = x > maxx ? x : maxx;

                miny = y < miny ? y : miny;
                maxy = y > maxy ? y : maxy;

                minz = z < minz ? z : minz;
                maxz = z > maxz ? z : maxz;
            }

            mx /= mtot;
            my /= mtot;
            mz /= mtot;

            maxx -= minx;
            maxy -= miny;
            maxz -= minz;

            int axis_to_split = 0;
            split_pos = mx;

            if (maxy > maxx)
            {
                axis_to_split = 1;
                split_pos = my;
                if (maxz > maxy)
                {
                    axis_to_split = 2;
                    split_pos = mz;
                }
            }
            else if (maxz > maxx)
            {
                axis_to_split = 2;
                split_pos = mz;
            }
            return axis_to_split;
        }

        Biter splitBodiesAlong(Biter begin, Biter end, unsigned int i, double split_pos)
        {
            if (std::distance(begin, end) > para_thres)
            {
                return __gnu_parallel::partition(begin, end,
                                                 std::bind([](const Body &b, int d, double pivot)
                                                           { return b.X[d] < pivot; },
                                                           std::placeholders::_1,
                                                           i,
                                                           split_pos));
            }
            else
            {
                return std::partition(begin, end,
                                      std::bind([](const Body &b, int d, double pivot)
                                                { return b.X[d] < pivot; },
                                                std::placeholders::_1,
                                                i,
                                                split_pos));
            }
        }

        void sortBodiesAlong(Biter begin, Biter end, unsigned int i)
        {
            auto mid = begin + std::distance(begin, end) / 2;

            if (std::distance(begin, end) > para_thres)
            {
                __gnu_parallel::nth_element(begin, mid, end,
                                            std::bind([](const Body &l, const Body &r, int d)
                                                      { return l.X[d] < r.X[d]; },
                                                      std::placeholders::_1,
                                                      std::placeholders::_2,
                                                      i));
            }
            else

            {
                std::nth_element(begin, mid, end,
                                 std::bind([](const Body &l, const Body &r, int d)
                                           { return l.X[d] < r.X[d]; },
                                           std::placeholders::_1,
                                           std::placeholders::_2,
                                           i));
            }
        }

        void buildCells(Biter begin, Biter end, size_t index, int level = 0)
        {
            auto &cell = cells[index];
#ifdef DEBUG
            std::cout << "build cell " << index << std::endl;
#endif
            cell.index = index;
            cell.BODY = std::distance(bodies.begin(), begin);
            cell.NBODY = std::distance(begin, end);
#ifdef DEBUG
            std::cout << cell.BODY << " " << cell.NBODY << "\n";
#endif
            cell.left = -1, cell.right = -1;

            if (cell.NBODY < ncrit)
            {
                return;
            }
            else
            {
                auto mid = begin + cell.NBODY / 2;
                int axis_to_split;
                double split_pos;

                if (treetype == kdtree)
                {
                    axis_to_split = (level + 1) % dim;
                    sortBodiesAlong(begin, end, axis_to_split);
                }
                else if (treetype == rcbtree)
                {
                    axis_to_split = findSplitAxis(begin, end, split_pos);
                    mid = splitBodiesAlong(begin, end, axis_to_split, split_pos);
                }

                cells.resize(cells.size() + 2);

                cell.left = static_cast<unsigned int>(cells.size() - 2);

                cell.right = static_cast<unsigned int>(cells.size() - 1);

                buildCells(begin, mid, cells[index].left, level + 1);

                buildCells(mid, end, cells[index].right, level + 1);
            }
        }

        void allocateMultipoles()
        {
            std::vector<std::vector<real_t>> Multipoles((int)cells.size(),
                                                        std::vector<real_t>(NTERM, 0));

            multipoles = std::move(Multipoles);

            std::vector<std::vector<real_t>> Locals((int)cells.size(),
                                                    std::vector<real_t>(NTERM, 0));

            locals = std::move(Locals);

            std::vector<std::vector<real_t>> Pns((int)cells.size(),
                                                 std::vector<real_t>((EXPANSION + 1), 0));
            pns = std::move(Pns);

            for (size_t i = 0; i < cells.size(); i++)
            {
                cells[i].M = &multipoles[i][0];
                cells[i].L = &locals[i][0];
                cells[i].Pn = &pns[i][0];
            }
        }

        void buildTree()
        {
            cells.resize(1);
            cells.reserve(bodies.size());

            buildCells(bodies.begin(), bodies.end(), 0, -1);
#ifdef DEBUG
            std::cout << "Tree construction done \n";
#endif
            allocateMultipoles();
#ifdef DEBUG
            std::cout << "Multipole and local memory allocation done \n";
#endif
        }

        void setBodies(size_t num)
        {
            Bodies bodies_new(num);
            bodies = std::move(bodies_new);
        }

        void preorder(size_t index)
        {
            auto &node = cells[index];
            count++;

            if (!node.isLeaf())
            {
                preorder(node.left);
                preorder(node.right);
            }
            else
            {
                P2M(&node);
                std::cout << node.X[0] << " " << node.X[1] << " " << node.X[2] << " " << node.R << " " << node.left << std::endl;
                return;
            }
        }

        void preOrder()
        {
            preorder(0);
        }

        void upwardPass(Cell *Ci)
        {
            if (Ci->isLeaf())
            {
                P2M(Ci);
            }
            else
            {
                for (Cell *Cj = &cells[Ci->left]; Cj != &cells[Ci->left] + 2; Cj++)
                {
#pragma omp task untied
                    upwardPass(Cj);
                }
#pragma omp taskwait

                M2M(Ci);
            }
        }

        //! Upward pass interface
        void upwardPass()
        {
#pragma omp parallel
#pragma omp single nowait

            upwardPass(&cells[0]);
        }

        // Set the tree to be either (0) KD tree (1) Recursive Coordinates Bisection tree
        void setType(int i)
        {
            if (i == 0)
            {
                treetype = kdtree;
            }
            else if (i == 1)
            {
                treetype = rcbtree;
            }
        }

        void P2M(Cell *C)
        {
            if (C->NBODY > 0)
            {
                Pointp P[C->NBODY];

                int cnt = 0;
                for (auto b = C->BODY; b != C->BODY + C->NBODY; b++, cnt++)
                {
                    auto B = bodies[b];
                    P[cnt] = Pointp(B.X[0], B.X[1], B.X[2]);
                }

                Min_spherep mb(P, P + C->NBODY);

                const double *center = mb.center_cartesian_begin();

                // use rminball as R
                C->X[0] = center[0];
                C->X[1] = center[1];
                C->X[2] = center[2];
                C->R = mb.radius();
            }
            else
            {
                C->R = 1e-6 * C->R;
            }

            real_t r_multipole[NTERM];

            for (int indice = 0; indice < NTERM; indice++)
                r_multipole[indice] = 0.0;

            for (auto b = C->BODY; b != C->BODY + C->NBODY; b++)
            {
                Body *B = &bodies[b];

                for (int d = 0; d < 3; d++)
                {
                    dX[d] = B->X[d] - C->X[d];
                }

                real_t Gnm[NTERM];

                make_Gnm_real(&dX[0], &Gnm[0], P);

                for (int indice = 0; indice < NTERM; indice++)
                    r_multipole[indice] += B->m * Gnm[indice];
            }

            for (int indice = 0; indice < NTERM; indice++)
                C->M[indice] += r_multipole[indice];
        }

        void M2M(Cell *Ci)
        {
            std::vector<Sphere> S;

            for (Cell *ci = &cells[Ci->left]; ci != &cells[Ci->left] + 2; ci++)
            {
                if (ci->left == -1)
                {
                    Point p(ci->X[0], ci->X[1], ci->X[2]);
                    S.push_back(Sphere(p, ci->R));
                }
                else
                {
                    // use granddaughers information to get R_max
                    for (Cell *cii = &cells[ci->left]; cii != &cells[ci->left] + 2;
                         cii++)
                    {
                        Point p(cii->X[0], cii->X[1], cii->X[2]);
                        S.push_back(Sphere(p, cii->R));
                    }
                }
            }

            Min_sphere mb(S.begin(), S.end());

            double rminball = mb.radius();

            const double *center = mb.center_cartesian_begin();
            Ci->X[0] = center[0];
            Ci->X[1] = center[1];
            Ci->X[2] = center[2];
            Ci->R = rminball;

            complex_t c_multipole[NTERM];
            for (int i = 0; i < NTERM; i++)
            {
                c_multipole[i] = 0;
            }

            for (Cell *Cj = &cells[Ci->left]; Cj != &cells[Ci->left] + 2; Cj++)
            {
                for (int d = 0; d < 3; d++)
                {
                    dX[d] = Cj->X[d] - Ci->X[d];
                }

                real_t r_multipole_Cj[NTERM];
                complex_t c_multipole_Cj[NTERM];

                for (int i = 0; i < NTERM; i++)
                {
                    r_multipole_Cj[i] = Cj->M[i];
                    c_multipole_Cj[i] = 0;
                }

                real_2_complex(r_multipole_Cj, c_multipole_Cj, P);

                complex_t Gnm[NTERM];

                make_Gnm(&dX[0], &Gnm[0], P);

                for (int n = 0; n <= P; n++)
                {
                    for (int m = 0; m <= n; m++)
                    {
                        for (int k = 0; k <= n; k++)
                            for (int l = std::max(-k, m - n + k); l <= std::min(k, m + n - k); l++)
                                c_multipole[index(n, m)] +=
                                    c_multipole_Cj[index(n - k, m - l)] * Gnm[index(k, l)];
                    }
                }
            }

            real_t r_multipole[NTERM];
            complex_2_real(c_multipole, r_multipole, P);

            for (int indice = 0; indice < NTERM; indice++)
                Ci->M[indice] += r_multipole[indice];
        }

        void M2L_rotate(Cell *Ci, Cell *Cj)
        {
        }

        void L2L(Cell *Ci)
        {
        }

        void L2P(Cell *Ci)
        {
        }

        void printTree()
        {
            for (auto &c : cells)
            {
                std::cout << c.index << " " << c.NBODY << " "
                          << c.X[0] << " " << c.X[1] << " " << c.X[2] << " " << c.R << " " << c.left << std::endl;
            }

            std::cout << "Multipoles of root cell: \n";
            for (int i = 0; i < NTERM; i++)
            {
                std::cout << cells[0].M[i] << " ";
            }
            std::cout << "\n";
        }

    private:
        Bodies bodies;
        Cells cells;

        std::vector<std::vector<real_t>> multipoles;
        std::vector<std::vector<real_t>> locals;
        std::vector<std::vector<real_t>> pns;

        TreeType treetype = kdtree;
        int dim = 3;
        int count = 0;
    };
}
