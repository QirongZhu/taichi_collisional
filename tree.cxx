#include "tree.h"

namespace FMM
{
  void Tree::initKernel()
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

  int Tree::findSplitAxis(Biter begin, Biter end, double &split_pos)
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

  Biter Tree::splitBodiesAlong(Biter begin, Biter end, unsigned int i, double split_pos)
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

  void Tree::sortBodiesAlong(Biter begin, Biter end, unsigned int i)
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

  void Tree::buildCells(Biter begin, Biter end, size_t index, int level = 0)
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
    cell.left = -1, cell.NCHILD = 0;

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

      cell.NCHILD = 2; //static_cast<unsigned int>(cells.size() - 1);

      buildCells(begin, mid, cells[index].left, level + 1);

      buildCells(mid, end, cells[index].left + 1, level + 1);
    }
  }

  void Tree::allocateMultipoles()
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

  void Tree::buildTree()
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

  void Tree::setBodies(Bodies &bodies_)
  {
    bodies = std::move(bodies_);
  }

  void Tree::preorder(size_t index)
  {
    auto &node = cells[index];
    count++;

    if (!node.isLeaf())
    {
      preorder(node.left);
      preorder(node.left+1);
    }
    else
    {
      P2M(&node);
      std::cout << node.X[0] << " " << node.X[1] << " " << node.X[2] << " " << node.R << " " << node.left << std::endl;
      return;
    }
  }

  void Tree::preOrder()
  {
    preorder(0);
  }

  // Set the tree to be either (0) KD tree (1) Recursive Coordinates Bisection tree
  void Tree::setType(int i)
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

  void Tree::printTree()
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

}
