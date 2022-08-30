#include "tree.h"

#include <vectorclass/vectorclass.h>

#ifdef DOUBLE_P2P
const int NSIMD = 8;
#else
const int NSIMD = 16;
#endif

namespace FMM
{
    inline int oddOrEven(int n)
    {
        return (((n)&1) == 1) ? -1 : 1;
    }

    real_t norm(real_t *X) { return X[0] * X[0] + X[1] * X[1] + X[2] * X[2]; }

    int index(int n, int m) { return n * (n + 1) + m; }

    void real_2_complex(real_t *real_arr, complex_t *complex_arr, int order)
    {
        for (int n = 0; n <= order; n++)
        {
            int index_start = n * n + n;
            complex_arr[index_start] = complex_t(real_arr[index_start], 0);
            real_t oddevenfac = -1;
            for (int m = 1; m <= n; m++)
            {
                complex_arr[index_start - m] = complex_t(real_arr[index_start + m] * oddevenfac,
                                                         -real_arr[index_start - m] * oddevenfac);
                complex_arr[index_start + m] = complex_t(real_arr[index_start + m], real_arr[index_start - m]);
                oddevenfac *= -1;
            }
        }
    }

    void complex_2_real(complex_t *complex_arr, real_t *real_arr, int order)
    {
        for (int n = 0; n <= order; n++)
        {
            int index_start = n * n + n;
            real_arr[index_start] = std::real(complex_arr[index_start]);
            for (int m = 1; m <= n; m++)
            {
                real_arr[index_start - m] = std::imag(complex_arr[index_start + m]);
                real_arr[index_start + m] = std::real(complex_arr[index_start + m]);
            }
        }
    }

    void make_Tnm(real_t *dX, complex_t *Tnm, int Order)
    {
        real_t r2 = norm(dX);

        real_t invr2 = 1 / r2;

        Tnm[0] = complex_t(sqrt(invr2), 0);

        for (int n = 1; n <= Order; n++)
        {
            int n2 = n * n;
            int ind = n2 + n + n;

            Tnm[ind] = Tnm[n2 - 1] * complex_t(dX[0] * (2 * n - 1) * invr2, dX[1] * (2 * n - 1) * invr2);

            if (n > 0)
            {
                int m, indexstart = n2 - n, indexnew = n2 + n;

                if (n > 1)
                {
                    for (m = 0; m < n - 1; m++)
                    {
                        Tnm[indexnew + m] =
                            invr2 * (2 * n - 1.) * dX[2] * Tnm[indexstart + m] -
                            invr2 * ((n - 1) * (n - 1) * (n > 1) - m * m) * Tnm[n2 - 3 * n + 2 + m];
                    }
                }

                m = n - 1;
                Tnm[indexnew + m] = invr2 * (2 * n - 1.) * dX[2] * Tnm[indexstart + m];
            }

            real_t oddoreven = -1;

            for (int m = 1; m <= n; m++)
            {
                Tnm[n2 + n - m] = std::conj(Tnm[n2 + n + m] * oddoreven);
                oddoreven *= -1;
            }
        }
    }

    void make_Gnm(real_t *dX, complex_t *Gnm, int Order)
    {
        real_t r2 = norm(dX);

        Gnm[0] = complex_t(1, 0);

        for (int n = 1; n <= Order; n++)
        {
            int n2 = n * n;
            int ind = n2 + n + n;

            Gnm[ind] = Gnm[n2 - 1] * complex_t(dX[0] / (2 * n), dX[1] / (2 * n));

            int m, m2, indexstart = n2 - n, indexnew = n2 + n;

            for (m = 0; m < n - 1; m++)
            {
                m2 = n2 - m * m;
                Gnm[indexnew + m] =
                    (2 * n - 1.) / m2 * dX[2] * Gnm[indexstart + m] - r2 / m2 * Gnm[n2 - 3 * n + 2 + m];
            }

            Gnm[indexnew + n - 1] = dX[2] * Gnm[indexstart + n - 1];

            real_t oddoreven = -1;
            for (int m = 1; m <= n; m++)
            {
                Gnm[indexnew - m] = std::conj(Gnm[indexnew + m] * oddoreven);
                oddoreven *= -1;
            }
        }
    }

    void make_Gnm_real(real_t *dX, real_t *Gnm, int Order)
    {
        real_t x = dX[0], y = dX[1], z = dX[2];
        real_t r2 = x * x + y * y + z * z;
        Gnm[0] = 1;
        for (int n = 1; n <= Order; n++)
        {
            Gnm[index(n, -(n - 1))] = z * Gnm[index(n - 1, -(n - 1))];
            Gnm[index(n, n - 1)] = z * Gnm[index(n - 1, n - 1)];

            for (int m = -(n - 2); m < (n - 1); m++)
                Gnm[index(n, m)] = ((2 * n - 1) * z * Gnm[index(n - 1, m)] -
                                    r2 * Gnm[index(n - 2, m)]) /
                                   (n * n - m * m);

            Gnm[index(n, -n)] = (Gnm[index(n - 1, n - 1)] * y + (n > 1) * Gnm[index(n - 1, -(n - 1))] * x) / (2.0 * n);
            Gnm[index(n, n)] = (Gnm[index(n - 1, n - 1)] * x - (n > 1) * Gnm[index(n - 1, -(n - 1))] * y) / (2.0 * n);
        }
    }

    void swap_x_z(real_t *arr)
    {
        real_t arrold[NTERM];
        memcpy(arrold, arr, NTERM * sizeof(real_t));

        arr[2] = arrold[3];
        arr[3] = arrold[2];

        real_t coeff;

#if EXPANSION > 1
        coeff = 1 / 2.0;
        arr[4] = (4 * arrold[5]) * coeff;
        arr[5] = (arrold[4]) * coeff;
        arr[6] = (-arrold[6] + arrold[8]) * coeff;
        /*arr[7] = (2*arrold[7])*coeff; */
        arr[8] = (3 * arrold[6] + arrold[8]) * coeff;
#endif

#include "multipole_rotation_matrix/rotate_order_3.txt"
#include "multipole_rotation_matrix/rotate_order_4.txt"
#include "multipole_rotation_matrix/rotate_order_5.txt"
#include "multipole_rotation_matrix/rotate_order_6.txt"
#include "multipole_rotation_matrix/rotate_order_7.txt"
#include "multipole_rotation_matrix/rotate_order_8.txt"
#include "multipole_rotation_matrix/rotate_order_9.txt"
#include "multipole_rotation_matrix/rotate_order_10.txt"
#include "multipole_rotation_matrix/rotate_order_11.txt"
#include "multipole_rotation_matrix/rotate_order_12.txt"
#include "multipole_rotation_matrix/rotate_order_13.txt"
#include "multipole_rotation_matrix/rotate_order_14.txt"
#include "multipole_rotation_matrix/rotate_order_15.txt"
#include "multipole_rotation_matrix/rotate_order_16.txt"
#include "multipole_rotation_matrix/rotate_order_17.txt"
#include "multipole_rotation_matrix/rotate_order_18.txt"
#include "multipole_rotation_matrix/rotate_order_19.txt"
#include "multipole_rotation_matrix/rotate_order_20.txt"
#include "multipole_rotation_matrix/rotate_order_21.txt"
#include "multipole_rotation_matrix/rotate_order_22.txt"
#include "multipole_rotation_matrix/rotate_order_23.txt"
#include "multipole_rotation_matrix/rotate_order_24.txt"
#include "multipole_rotation_matrix/rotate_order_25.txt"
#include "multipole_rotation_matrix/rotate_order_26.txt"
#include "multipole_rotation_matrix/rotate_order_27.txt"
#include "multipole_rotation_matrix/rotate_order_28.txt"
#include "multipole_rotation_matrix/rotate_order_29.txt"
#include "multipole_rotation_matrix/rotate_order_30.txt"
    }

    void Tree::P2M(Cell *C)
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

    void Tree::M2M(Cell *Ci)
    {
        std::vector<Sphere> S;

        for (Cell *ci = &cells[Ci->CHILD]; ci != &cells[Ci->CHILD] + Ci->NCHILD; ci++)
        {
            if (ci->CHILD == -1)
            {
                Point p(ci->X[0], ci->X[1], ci->X[2]);
                S.push_back(Sphere(p, ci->R));
            }
            else
            {
                // use granddaughers information to get R_max
                for (Cell *cii = &cells[ci->CHILD]; cii != &cells[ci->CHILD] + Ci->NCHILD;
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

        for (Cell *Cj = &cells[Ci->CHILD]; Cj != &cells[Ci->CHILD] + Ci->NCHILD; Cj++)
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

    void Tree::P2P(Cell *Ci, Cell *Cj)
    {
        Body *Bi = &bodies[Ci->BODY];
        Body *Bj = &bodies[Cj->BODY];

        int ni = Ci->NBODY;
        int nj = Cj->NBODY;

#if SIMD_P2P
        if (ni > NSIMD / 2)
        {
            int nii = (ni + NSIMD - 1) & (-NSIMD);

#ifndef DOUBLE_P2P
            float Xi[nii], Yi[nii], Zi[nii], Mi[nii];
            // float Xi[nii] __attribute__((aligned(64)));
            // float Yi[nii] __attribute__((aligned(64)));
            // float Zi[nii] __attribute__((aligned(64)));
            // float Mi[nii] __attribute__((aligned(64)));
#else
            double Xi[nii], Yi[nii], Zi[nii], Mi[nii];
            // double Xi[nii] __attribute__((aligned(64)));
            // double Yi[nii] __attribute__((aligned(64)));
            // double Zi[nii] __attribute__((aligned(64)));
            // double Mi[nii] __attribute__((aligned(64)));
#endif

            for (int k = 0; k < ni; k++)
            {
                Xi[k] = -Bi[k].X[0];
                Yi[k] = -Bi[k].X[1];
                Zi[k] = -Bi[k].X[2];
            }

#ifndef DOUBLE_P2P
            Vec16f xi, yi, zi, r, mj, mi;
            Vec16f invR, dx, dy, dz, r2;
            Vec16f ax, ay, az, pot;
#else
            Vec8d xi, yi, zi, r, mj, mi;
            Vec8d invR, dx, dy, dz, r2;
            Vec8d ax, ay, az, pot;
#endif
            for (int i = 0; i < nii; i = i + NSIMD)
            {
                xi.load(Xi + i);
                yi.load(Yi + i);
                zi.load(Zi + i);

                ax = 0, ay = 0, az = 0, pot = 0;

                for (int j = 0; j < nj; j++)
                {

                    mj = Bj[j].m;
                    dx = Bj[j].X[0] + xi;
                    dy = Bj[j].X[1] + yi;
                    dz = Bj[j].X[2] + zi;
                    r2 = dx * dx + dy * dy + dz * dz;

                    r = sqrt(r2);
                    invR = 1 / r;
                    invR = select(r2 > 0, invR, 0);
                    mj *= invR;

                    pot += mj;

                    mj = mj * (invR * invR);
                    ax += dx * mj;
                    ay += dy * mj;
                    az += dz * mj;
                }

                pot.store(Mi + i);
                ax.store(Xi + i);
                ay.store(Yi + i);
                az.store(Zi + i);
            }

            for (int i = 0; i < ni; i++)
            {

#pragma omp atomic
                Bi[i].p += (real_t)Mi[i];
#pragma omp atomic
                Bi[i].F[0] += (real_t)Xi[i];
#pragma omp atomic
                Bi[i].F[1] += (real_t)Yi[i];
#pragma omp atomic
                Bi[i].F[2] += (real_t)Zi[i];
            }
        }
        else
        {
            // do not use vectorized version if n1*n2 too small
            for (int i = 0; i < ni; i++)
            {

                real_t acc[3] = {0.0, 0.0, 0.0};
                real_t dx[3] = {0.0, 0.0, 0.0};

                real_t pot = 0;

                for (int j = 0; j < nj; j++)
                {

                    for (int d = 0; d < 3; d++)
                        dx[d] = Bj[j].X[d] - Bi[i].X[d];

                    real_t R2 = (dx[0] * dx[0] + dx[1] * dx[1] + dx[2] * dx[2]);

                    if (R2 > 0)
                    {
                        real_t invR2 = 1 / R2;
                        real_t invR = Bj[j].m * sqrt(invR2);
                        pot += invR;
                        for (int d = 0; d < 3; d++)
                        {
                            dx[d] *= invR2 * invR;
                            acc[d] += dx[d];
                        }
                    }
                }
#pragma omp atomic
                Bi[i].p += (real_t)pot;
#pragma omp atomic
                Bi[i].F[0] += (real_t)acc[0];
#pragma omp atomic
                Bi[i].F[1] += (real_t)acc[1];
#pragma omp atomic
                Bi[i].F[2] += (real_t)acc[2];
            }
        }
#else
        {
            // do not use vectorized version if n1*n2 too small
            for (int i = 0; i < ni; i++)
            {

                real_t acc[3] = {0.0, 0.0, 0.0};
                real_t dx[3] = {0.0, 0.0, 0.0};

                real_t pot = 0;

                for (int j = 0; j < nj; j++)
                {

                    for (int d = 0; d < 3; d++)
                        dx[d] = Bj[j].X[d] - Bi[i].X[d];

                    real_t R2 = (dx[0] * dx[0] + dx[1] * dx[1] + dx[2] * dx[2]);

                    if (R2 > 0)
                    {
                        real_t invR2 = 1 / R2;
                        real_t invR = Bj[j].m * sqrt(invR2);
                        pot += invR;
                        for (int d = 0; d < 3; d++)
                        {
                            dx[d] *= invR2 * invR;
                            acc[d] += dx[d];
                        }
                    }
                }
#pragma omp atomic
                Bi[i].p += (real_t)pot;
#pragma omp atomic
                Bi[i].F[0] += (real_t)acc[0];
#pragma omp atomic
                Bi[i].F[1] += (real_t)acc[1];
#pragma omp atomic
                Bi[i].F[2] += (real_t)acc[2];
            }
        }
#endif
    }

    //! Recursive call to dual tree traversal for horizontal pass
    void Tree::horizontalPass(Cell *Ci, Cell *Cj)
    {
        for (int d = 0; d < 3; d++)
            dX[d] = Ci->X[d] - Cj->X[d]; // Distance vector from source to target

        real_t R2 = norm(dX) * theta * theta; // Scalar distance squared

        if (R2 > (Ci->R + Cj->R) * (Ci->R + Cj->R))
        {                       // If distance is far enough
            M2L_rotate(Ci, Cj); //  M2L kernel
        }
        else if (Ci->isLeaf() && Cj->isLeaf())
        {                // Else if both cells are leafs
            P2P(Ci, Cj); //  P2P kernel
        }
        else if (Cj->isLeaf() || (Ci->R >= Cj->R && !Ci->isLeaf()))
        { // If Cj is leaf or Ci is larger
            for (Cell *ci = &cells[Ci->CHILD]; ci != &cells[Ci->CHILD] + Ci->NCHILD; ci++)
            {                                 // Loop over Ci's children
#pragma omp task untied if (ci->NBODY > 1000) //   Start OpenMP task if large enough task
                horizontalPass(ci, Cj);       //   Recursive call to target child cells
            }                                 //  End loop over Ci's children
        }
        else
        { // Else if Ci is leaf or Cj is larger
            for (Cell *cj = &cells[Cj->CHILD]; cj != &cells[Cj->CHILD] + Cj->NCHILD; cj++)
            {                           // Loop over Cj's children
                horizontalPass(Ci, cj); //   Recursive call to source child cells
            }                           //  End loop over Cj's children
        }                               // End if for leafs and Ci Cj size
    }

    //! Horizontal pass interface
    void Tree::horizontalPass()
    {
#pragma omp parallel                          // Start OpenMP
#pragma omp single nowait                     // Start OpenMP single region with nowait
        horizontalPass(&cells[0], &cells[0]); // Pass root cell to recursive call
    }

    void Tree::P2M_low(Cell *C)
    {
        real_t comx = 0, comy = 0, comz = 0, totalq = 0;

        for (auto b = C->BODY; b != C->BODY + C->NBODY; b++)
        {
            Body *B = &bodies[b];
            comx += B->X[0] * B->m;
            comy += B->X[1] * B->m;
            comz += B->X[2] * B->m;
            totalq += B->m;
        }

        comx /= totalq;
        comy /= totalq;
        comz /= totalq;

        C->X[0] = comx;
        C->X[1] = comy;
        C->X[2] = comz;

        real_t max_r2 = 1e-12 * C->R * C->R;

        for (auto b = C->BODY; b != C->BODY + C->NBODY; b++)
        {
            Body *B = &bodies[b];

            for (int d = 0; d < 3; d++)
            {
                dX[d] = B->X[d] - C->X[d];
            }

            max_r2 = std::max(max_r2, norm(dX));
        }

        C->R = sqrt(max_r2);
        C->M[0] = totalq;
    }

    void Tree::M2M_low(Cell *Ci)
    {
        real_t comx = 0, comy = 0, comz = 0, totalq = 0;

        for (Cell *Cj = &cells[Ci->CHILD]; Cj != &cells[Ci->CHILD] + Ci->NCHILD; Cj++)
        {
            comx += Cj->X[0] * Cj->M[0];
            comy += Cj->X[1] * Cj->M[0];
            comz += Cj->X[2] * Cj->M[0];
            totalq += Cj->M[0];
        }

        comx /= totalq;
        comy /= totalq;
        comz /= totalq;

        Ci->X[0] = comx;
        Ci->X[1] = comy;
        Ci->X[2] = comz;
        Ci->M[0] = totalq;

        real_t max_r = 0;

        for (Cell *Cj = &cells[Ci->CHILD]; Cj != &cells[Ci->CHILD] + Ci->NCHILD; Cj++)
        {
            for (int d = 0; d < 3; d++)
            {
                dX[d] = Cj->X[d] - Ci->X[d];
            }

            max_r = std::max(max_r, sqrt(norm(dX)) + Cj->R);
        }

        Ci->R = max_r;
    }

    void Tree::M2L_low(Cell *Ci, Cell *Cj)
    {
        for (int d = 0; d < 3; d++)
            dX[d] = Ci->X[d] - Cj->X[d];
        real_t f = Cj->M[0] / norm(dX);
        Ci->L[0] += f;
    }

    void Tree::P2P_low(Cell *Ci, Cell *Cj)
    {
        Body *Bi = &bodies[Ci->BODY];
        Body *Bj = &bodies[Cj->BODY];

        int ni = Ci->NBODY;
        int nj = Cj->NBODY;

#if SIMD_P2P
        int nii = (ni + NSIMD - 1) & (-NSIMD);

#ifndef DOUBLE_P2P
        float Xi[nii] __attribute__((aligned(64)));
        float Yi[nii] __attribute__((aligned(64)));
        float Zi[nii] __attribute__((aligned(64)));
#else
        double Xi[nii] __attribute__((aligned(64)));
        double Yi[nii] __attribute__((aligned(64)));
        double Zi[nii] __attribute__((aligned(64)));
#endif

        for (int k = 0; k < ni; k++)
        {
            Xi[k] = Bi[k].X[0];
            Yi[k] = Bi[k].X[1];
            Zi[k] = Bi[k].X[2];
        }

#ifndef DOUBLE_P2P
        Vec16f xi, yi, zi, r, mj, factor1, dx, dy, dz, r2;
#else
        Vec8d xi, yi, zi, r, mj, factor1, dx, dy, dz, r2;
#endif

        for (int i = 0; i < nii; i = i + NSIMD)
        {
            xi.load(Xi + i);
            yi.load(Yi + i);
            zi.load(Zi + i);
            factor1 = 0;
            for (int j = 0; j < nj; j++)
            {
                mj = Bj[j].m;
                dx = xi - Bj[j].X[0];
                dy = yi - Bj[j].X[1];
                dz = zi - Bj[j].X[2];
                r2 = dx * dx + dy * dy + dz * dz;
                mj = select(r2 > 0, mj, 0);
                r2 = select(r2 > 0, r2, HUGE);
                factor1 += mj / r2;
            }

            for (int k = 0; (k < NSIMD) && (i + k < ni); k++)
            {
                Bi[i + k].acc_old += (real_t)factor1[k];
            }
        }
#else
        for (int i = 0; i < ni; i++)
        {
            real_t invR2 = 0;

            for (int j = 0; j < nj; j++)
            {
                for (int d = 0; d < 3; d++)
                    dX[d] = Bi[i].X[d] - Bj[j].X[d];

                real_t R2 = norm(dX);

                if (R2 > 0)
                    invR2 += Bj[j].m / R2;
            }

            Bi[i].acc_old += invR2;
        }
#endif
    }

    void Tree::L2L_low(Cell *Ci)
    {
        for (Cell *Cj = &cells[Ci->CHILD]; Cj != &cells[Ci->CHILD] + Ci->NCHILD; Cj++)
        {
            Cj->L[0] += Ci->L[0];
        }
    }

    void Tree::L2P_low(Cell *C)
    {
        for (auto b = C->BODY; b != C->BODY + C->NBODY; b++)
        {
            bodies[b].acc_old += C->L[0];
        }
    }

    void Tree::upwardPass(Cell *Ci)
    {
        if (Ci->isLeaf())
        {
            P2M(Ci);
        }
        else
        {
            for (Cell *Cj = &cells[Ci->CHILD]; Cj != &cells[Ci->CHILD] + Ci->NCHILD; Cj++)
            {
#pragma omp task untied
                upwardPass(Cj);
            }
#pragma omp taskwait

            M2M(Ci);
        }
    }

    //! Upward pass interface
    void Tree::upwardPass()
    {
#pragma omp parallel
#pragma omp single nowait

        upwardPass(&cells[0]);
    }

    void Tree::upwardPass_low(Cell *Ci)
    {
        if (Ci->isLeaf())
        {
            P2M_low(Ci);
        }
        else
        {
            for (Cell *Cj = &cells[Ci->CHILD]; Cj != &cells[Ci->CHILD] + Ci->NCHILD; Cj++)
            {
#pragma omp task untied
                upwardPass_low(Cj);
            }
#pragma omp taskwait

            M2M_low(Ci);
        }
    }

    //! Upward pass interface
    void Tree::upwardPass_low()
    {
#pragma omp parallel
#pragma omp single nowait

        upwardPass_low(&cells[0]);
    }

    void Tree::M2L_rotate(Cell *Ci, Cell *Cj)
    {
        // if(Ci->NBODY <= 8){
        //   M2P(Ci, Cj);
        //   return;
        // }

        // if(Cj->NBODY <= 8)
        //   {
        // P2L(Ci, Cj);
        // return;
        //   }

        for (int d = 0; d < 3; d++)
            dX[d] = Ci->X[d] - Cj->X[d];

        real_t r2_xy = dX[0] * dX[0] + dX[1] * dX[1];
        real_t r_xy = sqrt(r2_xy);
        real_t r2 = r2_xy + dX[2] * dX[2];
        real_t r = sqrt(r2);
        real_t r_inv = 1.0 / r;

        real_t exp_a_z_re[P + 1], exp_a_z_im[P + 1];
        real_t exp_a_x_re[P + 1], exp_a_x_im[P + 1];

        exp_a_z_re[0] = 1;
        exp_a_z_im[0] = 0;
        exp_a_x_re[0] = 1;
        exp_a_x_im[0] = 0;
        exp_a_z_re[1] = 1;
        exp_a_z_im[1] = 0;
        exp_a_x_re[1] = dX[2] * r_inv;
        exp_a_x_im[1] = -r_xy * r_inv;

        if (r_xy > 0)
        {
            r_xy = 1.0 / r_xy;
            exp_a_z_re[1] = dX[1] * r_xy;
            exp_a_z_im[1] = dX[0] * r_xy;
        }

        for (int m = 1; m < P; m++)
        {
            exp_a_z_re[m + 1] = exp_a_z_re[m] * exp_a_z_re[1] - exp_a_z_im[m] * exp_a_z_im[1];
            exp_a_z_im[m + 1] = exp_a_z_re[m] * exp_a_z_im[1] + exp_a_z_im[m] * exp_a_z_re[1];
            exp_a_x_re[m + 1] = exp_a_x_re[m] * exp_a_x_re[1] - exp_a_x_im[m] * exp_a_x_im[1];
            exp_a_x_im[m + 1] = exp_a_x_re[m] * exp_a_x_im[1] + exp_a_x_im[m] * exp_a_x_re[1];
        }

        real_t r_multipole_Cj[NTERM], r_local_ci[NTERM];

        // Step 0: cache multipoles, scale them by r_invs, and make them into homogenius polynomials

        real_t scale_fac = 1;
        for (int n = 0; n <= P; n++)
        {
            for (int m = -n; m <= n; m++)
            {
                r_local_ci[index(n, m)] = 0;
                r_multipole_Cj[index(n, m)] = Cj->M[index(n, m)] * factorial_coef_oned[index(n, m)] * scale_fac;
            }
            scale_fac *= r_inv;
        }

        // Step 1:       first rotate multipoles around the original z - axis with alpha_z
        for (int n = 1; n <= P; n++)
        {
            int index_start = n * n + n;
            real_t Re_ei, Im_ei, Im_r, Re_r;
            for (int m = 1; m <= n; m++)
            {
                Re_ei = exp_a_z_re[m];
                Im_ei = exp_a_z_im[m];
                Im_r = r_multipole_Cj[index_start - m];
                Re_r = r_multipole_Cj[index_start + m];
                r_multipole_Cj[index_start - m] = Re_ei * Im_r + Im_ei * Re_r;
                r_multipole_Cj[index_start + m] = Re_ei * Re_r - Im_ei * Im_r;
            }
        }

        // Step 2:       swap x and z
        swap_x_z(r_multipole_Cj);

        // step 3:       rotate around(new) z - axis by alpha_x
        for (int n = 1; n <= P; n++)
        {
            int index_start = n * n + n;
            real_t Re_ei, Im_ei, Im_r, Re_r;
            for (int m = 1; m <= n; m++)
            {
                Re_ei = exp_a_x_re[m];
                Im_ei = exp_a_x_im[m];
                Im_r = r_multipole_Cj[index_start - m];
                Re_r = r_multipole_Cj[index_start + m];
                r_multipole_Cj[index_start - m] = Re_ei * Im_r + Im_ei * Re_r;
                r_multipole_Cj[index_start + m] = Re_ei * Re_r - Im_ei * Im_r;
            }
        }

        // Step 4:       swap x and z
        swap_x_z(r_multipole_Cj);

        // Lastly multiply multipoles by proper normalizations
        for (int n = 0; n < NTERM; n++)
        {
            r_multipole_Cj[n] *= factorial_coef_inv_oned[n];
        }

        // Actual computation of M2L, using the fact that along (0, 0, 1) vector
        // the singular solid harmonics is much simpler than the full form

        for (int n = 0; n <= P; n++)
        {
            int index_start = n * (n + 1);

#ifndef DOUBLEHEIGHT
            for (int k = 0; k <= P - n; k++)
#else
            for (int k = 0; k <= P; k++)
#endif
            {
                real_t fac = factorial_table[n + k];
                r_local_ci[index_start] += r_multipole_Cj[k * (k + 1)] * fac;
            }

            real_t fac;

#ifndef DOUBLEHEIGHT
            for (int m = 1; m <= std::min(n, P - n); m++)
#else
            for (int m = 1; m <= n; m++)
#endif
            {

#ifndef DOUBLEHEIGHT
                for (int k = m; k <= P - n; k++)
#else
                for (int k = m; k <= P; k++)
#endif
                {
                    int index_start_k = index(k, 0);
                    fac = factorial_table[(n + k)] * oddOrEven(m);
                    r_local_ci[index_start - m] += r_multipole_Cj[index_start_k - m] * fac;
                    r_local_ci[index_start + m] += r_multipole_Cj[index_start_k + m] * fac;
                }
            }
        }

        // Step 5:       swap x and z
        swap_x_z(r_local_ci);

        // Step 6:       Rotate local expansion around z - axis by(-alpha_x)
        for (int n = 1; n <= P; n++)
        {
            int index_start = n * n + n;
            real_t Re_ei, Im_ei, Im_c, Re_c;
            for (int m = 1; m <= n; m++)
            {
                Re_ei = exp_a_x_re[m];
                Im_ei = exp_a_x_im[m];
                Im_c = r_local_ci[index_start - m];
                Re_c = r_local_ci[index_start + m];
                r_local_ci[index_start - m] = Re_ei * Im_c - Im_ei * Re_c;
                r_local_ci[index_start + m] = Re_ei * Re_c + Im_ei * Im_c;
            }
        }

        // Step 7:       swap x and z
        swap_x_z(r_local_ci);

        // Step 8:       Final rotation around z - axis by(-alpha_z)
        for (int n = 1; n <= P; n++)
        {
            int index_start = n * n + n;
            real_t Re_ei, Im_ei, Im_c, Re_c;
            for (int m = 1; m <= n; m++)
            {
                Re_ei = exp_a_z_re[m];
                Im_ei = exp_a_z_im[m];
                Im_c = r_local_ci[index_start - m];
                Re_c = r_local_ci[index_start + m];
                r_local_ci[index_start - m] = Re_ei * Im_c - Im_ei * Re_c;
                r_local_ci[index_start + m] = Re_ei * Re_c + Im_ei * Im_c;
            }
        }

        // Step 9: scale Lmn by r_inv factors
        scale_fac = r_inv;
        for (int n = 0; n <= P; n++)
        {
            for (int m = -n; m <= n; m++)
            {
                r_local_ci[index(n, m)] *= scale_fac;
            }
            scale_fac *= r_inv;
        }

        // Now local expansion is expressed consistently with the original coordinates
        //  Store the results in real - valued form
        for (int indice = 0; indice < NTERM; indice++)
        {
#pragma omp atomic
            Ci->L[indice] += r_local_ci[indice];
        }
    }

}
