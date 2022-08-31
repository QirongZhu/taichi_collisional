#include "radixtree.h"
#include <bitset>

namespace FMM
{

    void RadixTree::setBodies(Bodies &bodies_)
    {
        bodies = std::move(bodies_);
    }

    int RadixTree::delta(int x, int y)
    {
        if (x >= 0 && x < sortedMortonCodes.size() && y >= 0 && y < sortedMortonCodes.size())
        {
            return CLZ(sortedMortonCodes[x].first ^ sortedMortonCodes[y].first);
        }
        return -1;
    }

    int2 RadixTree::determineRange(int idx)
    {
        int d = sign(delta(idx, idx + 1) - delta(idx, idx - 1));

        int dmin = delta(idx, idx - d);

        int lmax = 2;
        while (delta(idx, idx + lmax * d) > dmin)
            lmax = lmax * 2;

        int l = 0;
        for (int t = lmax / 2; t >= 1; t /= 2)
        {
            if (delta(idx, idx + (l + t) * d) > dmin)
                l += t;
        }

        int j = idx + l * d;

        int lr, rr;
        lr = idx < j ? idx : j;
        rr = idx > j ? idx : j;

        int2 range = std::make_pair(lr, rr);
        return range;
    }

    int RadixTree::findSplit(int2 range)
    {
        int first = range.first;
        int last = range.second;

        HashType firstCode = sortedMortonCodes[first].first;
        HashType lastCode = sortedMortonCodes[last].first;

        if (firstCode == lastCode)
            return (first + last) >> 1;

        int commonPrefix = CLZ(firstCode ^ lastCode);

        // Use binary search to find where the next bit differs.
        // Specifically, we are looking for the highest object that
        // shares more than commonPrefix bits with the first one.

        int split = first; // initial guess
        int step = last - first;

        do
        {
            step = (step + 1) >> 1;      // exponential decrease
            int newSplit = split + step; // proposed new position

            if (newSplit < last)
            {
                HashType splitCode = sortedMortonCodes[newSplit].first;

                int splitPrefix = CLZ(firstCode ^ splitCode);

                if (splitPrefix > commonPrefix)
                    split = newSplit; // accept proposal
            }
        } while (step > 1);

        return split;
    }

    void RadixTree::assignNode(int idx)
    {
        int2 range = determineRange(idx);

        int split = findSplit(range);

        Node *ChA;
        Node *ChB;

        if (split == range.first)
            ChA = &Leafs[split];
        else
            ChA = &Tree[split];

        if (split + 1 == range.second)
            ChB = &Leafs[split + 1];
        else
            ChB = &Tree[split + 1];

        Tree[idx].left = ChA;
        Tree[idx].right = ChB;
    }

    void RadixTree::getBoundBox()
    {

        double xmax = -HUGE, ymax = -HUGE, zmax = -HUGE;
        double xmin = HUGE, ymin = HUGE, zmin = HUGE;
        for (auto &p : bodies)
        {
            if (p.X[0] > xmax)
                xmax = p.X[0];
            else if (p.X[0] < xmin)
                xmin = p.X[0];

            if (p.X[1] > ymax)
                ymax = p.X[1];
            else if (p.X[1] < ymin)
                ymin = p.X[1];

            if (p.X[2] > zmax)
                zmax = p.X[2];
            else if (p.X[2] < zmin)
                zmin = p.X[2];
        }
        range_x = std::make_pair(xmin, xmax);
        range_y = std::make_pair(ymin, ymax);
        range_z = std::make_pair(zmin, zmax);

        double radius = range_x.second - range_x.first;

        if (radius < range_y.second - range_y.first)
            radius = range_y.second - range_y.first;

        if (radius < range_z.second - range_z.first)
            radius = range_z.second - range_z.first;

        radius = radius * 1.001;
        R0 = radius / 2;

        X0[0] = (range_x.first + range_x.second) / 2 - R0;
        X0[1] = (range_y.first + range_y.second) / 2 - R0;
        X0[2] = (range_z.first + range_z.second) / 2 - R0;

        // std::cout << R0 << " " << X0[0] << " " << X0[1] << " " << X0[2] << std::endl;
    }

    void RadixTree::sortMontonID()
    {
        getBoundBox();

        sortedMortonCodes.resize(bodies.size());

        const int max_int = (1 << 20);

#pragma omp parallel for
        for (size_t b = 0; b < bodies.size(); b++)
        {
            double x_f = (bodies[b].X[0] - X0[0]) / (2.0 * R0);
            HashType x = (HashType)floor(x_f * (double)max_int);

            double y_f = (bodies[b].X[1] - X0[1]) / (2.0 * R0);
            HashType y = (HashType)floor(y_f * (double)max_int);

            double z_f = (bodies[b].X[2] - X0[2]) / (2.0 * R0);
            HashType z = (HashType)floor(z_f * (double)max_int);

            HashType xx = expandBits((HashType)x);
            HashType yy = expandBits((HashType)y);
            HashType zz = expandBits((HashType)z);

            sortedMortonCodes[b] = std::make_pair(xx * 4 + yy * 2 + zz, b);
        }

        __gnu_parallel::sort(begin(sortedMortonCodes), end(sortedMortonCodes),
                             [](const int2 &l, const int2 &r)
                             { return l.first < r.first; });

        Bodies tmp = bodies;

#pragma omp paralle for
        for (size_t b = 0; b < bodies.size(); b++)
        {
            bodies[b] = tmp[sortedMortonCodes[b].second];

            // std::bitset<32> mid(sortedMortonCodes[b].first);
            // std::cout << mid << " " << sortedMortonCodes[b].second;
            //  std::cout << bodies[b].index << " " << bodies[b].X[0];
            // std::cout << bodies[b].index << " " << bodies[b].X[0] << " "
            //           << bodies[b].X[1] << " " << bodies[b].X[2] << std::endl;
        }
    }

    void RadixTree::buildRadixTree()
    {
        sortMontonID();

        Leafs.resize(bodies.size());
        for (size_t b = 0; b < Leafs.size(); b++)
        {
            Leafs[b].idx = b;
            Leafs[b].isleaf = true;
        }

        Tree.resize(bodies.size() - 1);
#pragma omp parallel for
        for (size_t b = 0; b < Tree.size(); b++)
        {
            Tree[b].idx = b;
            Tree[b].isleaf = false;
            assignNode(b);
        }

        root = &Tree[0];
    }

    void RadixTree::traverse(Node *n)
    {

        std::cout << "node[" << n->idx << "] ->cell[" << start << "] ";
        std::cout << " start: " << n->BODY << " end: " << n->BODY + n->NBODY << " cnt:" << n->NBODY << " \n";

        n->index = start;

        if (n->NBODY >= ncrit)
        {
            ++start;
            traverse(n->left);

            ++start;
            traverse(n->right);
        }
    }

    void RadixTree::traverse()
    {
        traverse(&Tree[0]);
    }

    void RadixTree::upwardPass()
    {
        upwardPass(&Tree[0]);
    }

    void RadixTree::upwardPass(Node *n)
    {
        //   std::cout << "On node " << n->idx << " " << n->index << "\n";

        if (n->isLeaf())
        {
            n->BODY = n->idx;
            n->NBODY = 1;
        }
        else
        {
            if (n->left)
                upwardPass(n->left);

            if (n->right)
                upwardPass(n->right);

            n->BODY = (n->left->BODY < n->right->BODY) ? n->left->BODY : n->right->BODY;
            n->NBODY = n->left->NBODY + n->right->NBODY;
        }
    }

    void RadixTree::printTree()
    {
        for (size_t c = 0; c < Tree.size(); c++)
        {
            std::cout << Tree[c].idx << " " << std::endl;
        }
    }

}