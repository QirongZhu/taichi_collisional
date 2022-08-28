#pragma once

#include <iostream>
#include <vector>
#include <algorithm>
#include <functional>
#include <complex>
#include <random>
#include <string>
#include <chrono>
#include <numeric>
#include <parallel/algorithm>

namespace FMM
{
    typedef double real_t;                  //!< Real type
    typedef std::complex<real_t> complex_t; //!< Complex type

    std::random_device rd;
    std::mt19937 e2(rd());
    std::uniform_real_distribution<> dist(-10, 10);

    enum TreeType
    {
        kdtree,
        rcbtree
    };

    struct Body
    {
        real_t X[3];
        real_t m;
        real_t F[3];
        real_t p;
        size_t index;

        Body()
        {
            X[0] = dist(e2);
            X[1] = dist(e2);
            X[2] = dist(e2);
            m = 1.0;
        }
    };

    typedef std::vector<Body> Bodies;
    typedef Bodies::iterator Biter;

    struct Cell
    {
        real_t X[3];
        real_t R;

        real_t *M;
        real_t *L;
        real_t *Pn;
        real_t min_acc;

        int left;
        int right;
        size_t BODY;
        size_t NBODY;
        size_t index;

        bool has_sink;
        bool has_source;

        bool isLeaf() { return left == -1; }
    };

    typedef std::vector<Cell> Cells;

    const int P = EXPANSION;                             //!< Order of expansions
    const int NTERM = (EXPANSION + 1) * (EXPANSION + 1); //!< Number of coefficients
    const int ncrit = 192;                               //!< Number of bodies per leaf cell
    const real_t theta = 0.65;                           //!< Multipole acceptance criterion
    const int para_thres = 1500;

    real_t dX[3];             //!< Distance vector
#pragma omp threadprivate(dX) //!< Make global variables private

    static real_t factorial_table[2 * EXPANSION + 8];

    static real_t factorial_coef[EXPANSION + 1][EXPANSION + 1];
    static real_t factorial_coef_inv[EXPANSION + 1][EXPANSION + 1];

    static real_t factorial_coef_oned[(EXPANSION + 1) * (EXPANSION + 1)];
    static real_t factorial_coef_inv_oned[(EXPANSION + 1) * (EXPANSION + 1)];

    static real_t combinator_coef[EXPANSION + 1][EXPANSION + 1];

}