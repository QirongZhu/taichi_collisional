#pragma once

#include <iostream>
#include <vector>
#include <algorithm>
#include <functional>
#include <random>
#include <string>
#include <chrono>
#include <numeric>
#include <parallel/algorithm>

namespace FMM
{

    std::random_device rd;
    std::mt19937 e2(rd());
    std::uniform_real_distribution<> dist(-10, 10);

    enum TreeType
    {
        kdtree,
        rcbtree
    };

    const int thres = 256;
    const int para_thres = 1500;

    struct Body
    {
        double X[3];
        double m;

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

        double X[3];
        double R;

        int left;
        int right;
        size_t BODY;
        size_t NBODY;
        size_t index;

        bool isLeaf() { return left == -1; }

        // unsigned int level;
    };

    typedef std::vector<Cell> Cells;

}
