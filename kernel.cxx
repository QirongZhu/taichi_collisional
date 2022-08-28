#include "tree.h"

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


}