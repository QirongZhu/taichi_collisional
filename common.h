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
        real_t acc_old;
        Body(real_t x_ = 0, real_t y_ = 0, real_t z_ = 0, real_t m_ = 0, size_t id = 0)
        {
            X[0] = x_, X[1] = y_, X[2] = z_, m = m_, index = id;
            F[0] = 0, F[1] = 0, F[2] = 0, p = 0, acc_old = 0;
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
        int NCHILD;
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
    const int ncrit = 100;                               //!< Number of bodies per leaf cell
    const real_t theta = 0.4;                           //!< Multipole acceptance criterion
    const int para_thres = 1500;

    static real_t dX[3];      //!< Distance vector
#pragma omp threadprivate(dX) //!< Make global variables private

    static real_t factorial_table[2 * EXPANSION + 8];

    static real_t factorial_coef[EXPANSION + 1][EXPANSION + 1];
    static real_t factorial_coef_inv[EXPANSION + 1][EXPANSION + 1];

    static real_t factorial_coef_oned[NTERM];
    static real_t factorial_coef_inv_oned[NTERM];

    static real_t combinator_coef[EXPANSION + 1][EXPANSION + 1];

}