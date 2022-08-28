#include "tree.h"

#include <vectorclass/vectorclass.h>

const int NSIMD = 8;

namespace FMM
{
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

}