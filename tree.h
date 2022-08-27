#include <iostream>
#include <vector>
#include <algorithm>
#include <functional>
#include <random>
#include <string>
#include <chrono>

#include <numeric>

#include <parallel/algorithm>

std::random_device rd;
std::mt19937 e2(rd());
std::uniform_real_distribution<> dist(-10, 10);

enum TreeType{kdtree, rcbtree};

const int thres = 256;
const int para_thres = 5000;

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
  int left;
  int right;
  size_t part;
  size_t npart;
  size_t index;
  unsigned int level;
};

typedef std::vector<Cell> Cells;

class Tree
{
    
public:
    
    
  int findSplitAxis(Biter begin, Biter end, double& split_pos)
    {
      double mx=0, my=0, mz=0, mtot=0;
      double minx = 1e30, miny = 1e30, minz=1e30;
      double maxx = -1e30,maxy = -1e30, maxz = -1e30;
      auto it = begin;
      
      for(it =begin; it!=end; it++)
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
     
        return std::partition(begin, end,
                              std::bind([] (const Body&b, int d, double pivot)
                                        {return b.X[d] < pivot;},
                                    std::placeholders::_1,
                                    i,
                                    split_pos));
    }
    
    
  void sortBodiesAlong(Biter begin, Biter end, unsigned int i)
  {
    auto mid = begin + std::distance(begin, end)/2;
        
    if (std::distance(begin, end) > para_thres)
    {
        __gnu_parallel::nth_element(begin, mid, end,
		     std::bind([] (const Body&l, const Body&r, int d)
                       {return l.X[d] < r.X[d];},
			       std::placeholders::_1,
			       std::placeholders::_2,
			       i) );
    }
    else
     
      {
        std::nth_element(begin, mid, end,
                 std::bind([] (const Body&l, const Body&r, int d)
                           {return l.X[d] < r.X[d];},
                       std::placeholders::_1,
                       std::placeholders::_2,
                       i) );
      }
  }
    
  void buildCells(Biter begin, Biter end, size_t index, int level = 0)
  {
    auto & cell = cells[index];
#ifdef DEBUG
    std::cout << "build cell " << index << std::endl;
#endif
    cell.index = index;
    cell.part = std::distance(bodies.begin(), begin);
    cell.npart= std::distance(begin, end);
#ifdef DEBUG
    std::cout << cell.part << " " << cell.npart << "\n";
#endif
    cell.left = -1, cell.right = -1, cell.level = level;
        
    if (cell.npart < thres)
      {
	return;
      }
    else
      {
	auto mid = begin + cell.npart/2;
	int axis_to_split;
	double split_pos;
            
	if (treetype == kdtree)
	  {
	    axis_to_split = (level+1) % dim;
	    sortBodiesAlong(begin, end, axis_to_split);
	  }
    else if (treetype == rcbtree)
    {
        axis_to_split = findSplitAxis(begin, end, split_pos);
        mid = splitBodiesAlong(begin, end, axis_to_split, split_pos);
    }
            
	cells.resize(cells.size()+2);
                        
	cell.left = static_cast<unsigned int> (cells.size() - 2);
            
	cell.right = static_cast<unsigned int> (cells.size() - 1);
            
            
	buildCells(begin, mid, cells[index].left, level+1);
            
	buildCells(mid, end, cells[index].right, level+1);
            
      }
        
  }
    
  void buildTree()
  {
    cells.resize(1);
    cells.reserve(bodies.size());
    buildCells(bodies.begin(), bodies.end(), 0, -1);
  }
    
  void setBodies(size_t num)
  {
    Bodies bodies_new(num);
    bodies = std::move(bodies_new);
  }
    
  void preorder(size_t index)
  {
    auto& node = cells[index];
    count++;
        
    if (node.left != -1 && node.right != -1)
      {
	preorder(node.left);
	preorder(node.right);
      }
    else
      {
	return;
      }
        
  }
    
  void preOrder()
  {
    preorder(0);
  }
    
  void setType(int i)
  {
    if (i==0)
      {
	treetype = kdtree;
      }
    if (i==1)
      {
	treetype = rcbtree;
      }
  }
    
private:
  Bodies bodies;
  Cells cells;
  TreeType treetype = kdtree;
  int dim = 3;
  int count = 0;
};
