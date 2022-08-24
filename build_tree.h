#ifndef buildtree_h
#define buildtree_h
#include <algorithm>
#include "exafmm.h"
#include "timer.h"

#ifndef MAX_DEPTH
#define MAX_DEPTH 28
#endif

#include <omp.h>

namespace omp_par
{
  template < class T, class StrictWeakOrdering >
  void merge(T A_, T A_last, T B_, T B_last, T C_, size_t p, StrictWeakOrdering comp)
  {
    typedef typename std::iterator_traits < T >::difference_type _DiffType;
    typedef typename std::iterator_traits < T >::value_type _ValType;

    _DiffType N1 = A_last - A_;
    _DiffType N2 = B_last - B_;
        
    if(N1 == 0 && N2 == 0)
      return;
    if(N1 == 0 || N2 == 0)
      {
	_ValType *A = (N1 == 0 ? &B_[0] : &A_[0]);
	_DiffType N = (N1 == 0 ? N2 : N1);

#pragma omp parallel for
	for(size_t i = 0; i < p; i++)
	  {
	    _DiffType indx1 = (i * N) / p;
	    _DiffType indx2 = ((i + 1) * N) / p;
	    memcpy(&C_[indx1], &A[indx1], (indx2 - indx1) * sizeof(_ValType));
	  }
	return;
      }

    //Split both arrays ( A and B ) to n equal parts.
    //Find the position of each split in the final merged array.

    size_t n = 32;

    _ValType *split = new _ValType[p * n * 2];
    _DiffType *split_size = new _DiffType[p * n * 2];

#pragma omp parallel for
    for(size_t i = 0; i < p; i++)
      {
	for(size_t j = 0; j < n; j++)
	  {
	    size_t indx = i * n + j;
	    _DiffType indx1 = (indx * N1) / (p * n);

	    split[indx] = A_[indx1];
	    split_size[indx] = indx1 + (std::lower_bound(B_, B_last, split[indx], comp) - B_);

	    indx1 = (indx * N2) / (p * n);
	    indx += p * n;
	    split[indx] = B_[indx1];
	    split_size[indx] = indx1 + (std::lower_bound(A_, A_last, split[indx], comp) - A_);
	  }
      }

    //Find the closest split position for each thread that will
    //divide the final array equally between the threads.
    _DiffType *split_indx_A = new _DiffType[p + 1];
    _DiffType *split_indx_B = new _DiffType[p + 1];

    split_indx_A[0] = 0;
    split_indx_B[0] = 0;

    split_indx_A[p] = N1;
    split_indx_B[p] = N2;

#pragma omp parallel for
    for(size_t i = 1; i < p; i++)
      {
	_DiffType req_size = (i * (N1 + N2)) / p;

	size_t j = std::lower_bound(&split_size[0], &split_size[p * n], req_size,
				    std::less < _DiffType > ()) - &split_size[0];

	if(j >= p * n)
	  j = p * n - 1;
	_ValType split1 = split[j];
	_DiffType split_size1 = split_size[j];

	j =
	  (std::lower_bound(&split_size[p * n], &split_size[p * n * 2], req_size,
			    std::less < _DiffType > ()) - &split_size[p * n]) + p * n;
	if(j >= 2 * p * n)
	  j = 2 * p * n - 1;
	if(abs(split_size[j] - req_size) < abs(split_size1 - req_size))
	  {
	    split1 = split[j];
	    split_size1 = split_size[j];
	  }

	split_indx_A[i] = std::lower_bound(A_, A_last, split1, comp) - A_;
	split_indx_B[i] = std::lower_bound(B_, B_last, split1, comp) - B_;
      }
    delete[]split;
    delete[]split_size;

    //Merge for each thread independently.

#pragma omp parallel for
    for(size_t i = 0; i < p; i++)
      {
	T C = C_ + split_indx_A[i] + split_indx_B[i];

	std::merge(A_ + split_indx_A[i], A_ + split_indx_A[i + 1],
		   B_ + split_indx_B[i], B_ + split_indx_B[i + 1], C, comp);
      }
    delete[]split_indx_A;
    delete[]split_indx_B;
  }

  template < class T, class StrictWeakOrdering > void merge_sort(T A, T A_last, StrictWeakOrdering comp)
  {
    typedef typename std::iterator_traits < T >::difference_type _DiffType;
    typedef typename std::iterator_traits < T >::value_type _ValType;

    size_t p = omp_get_max_threads();

    _DiffType N = A_last - A;

    if(N < 2 * p)
      {
	std::sort(A, A_last, comp);
	return;
      }

    //Split the array A to p equal parts.
    _DiffType *split = new _DiffType[p + 1];

    split[p] = N;
#pragma omp parallel for
    for(size_t id = 0; id < p; id++)
      {
	split[id] = (id * N) / p;
      }

    //Sort each part independently.
#pragma omp parallel for
    for(size_t id = 0; id < p; id++)
      {
	std::sort(A + split[id], A + split[id + 1], comp);
      }

    //Merge two parts at a time.
    _ValType *B = new _ValType[N];
    _ValType *A_ = &A[0];
    _ValType *B_ = &B[0];

    for(size_t j = 1; j < p; j = j * 2)
      {
	for(size_t i = 0; i < p; i = i + 2 * j)
	  {
	    if(i + j < p)
	      {
		merge(A_ + split[i], A_ + split[i + j], A_ + split[i + j],
		      A_ + split[(i + 2 * j <= p ? i + 2 * j : p)], B_ + split[i], p, comp);
	      }
	    else
	      {
#pragma omp parallel for
		for(size_t k = split[i]; k < split[p]; k++)
		  B_[k] = A_[k];
	      }
	  }
	_ValType *tmp_swap = A_;

	A_ = B_;
	B_ = tmp_swap;
      }

    //The final result should be in A.
    if(A_ != &A[0])
      {
#pragma omp parallel for
	for(size_t i = 0; i < N; i++)
	  A[i] = A_[i];
      }

    //Free memory.
    delete[]split;
    delete[]B;
  }

  template < class T > void merge_sort(T A, T A_last)
  {
    typedef typename std::iterator_traits < T >::value_type _ValType;

    merge_sort(A, A_last, std::less < _ValType > ());
  }
}



namespace exafmm
{
  //! Get bounding box of bodies
  void getBounds(Bodies & bodies, real_t & R0, real_t * X0)
  {
    real_t Xmin[3], Xmax[3];	// Min, max of domain
    for(int d = 0; d < 3; d++)
      Xmin[d] = Xmax[d] = bodies[0].X[d];	// Initialize Xmin, Xmax
    for(size_t b = 0; b < bodies.size(); b++)
      {				// Loop over range of bodies
	for(int d = 0; d < 3; d++)
	  Xmin[d] = fmin(bodies[b].X[d], Xmin[d]);	//  Update Xmin
	for(int d = 0; d < 3; d++)
	  Xmax[d] = fmax(bodies[b].X[d], Xmax[d]);	//  Update Xmax
      }				// End loop over range of bodies
    for(int d = 0; d < 3; d++)
      X0[d] = (Xmax[d] + Xmin[d]) / 2;	// Calculate center of domain
    R0 = 0;			// Initialize localRadius
    for(int d = 0; d < 3; d++)
      {				// Loop over dimensions
	R0 = fmax(X0[d] - Xmin[d], R0);	//  Calculate min distance from center
	R0 = fmax(Xmax[d] - X0[d], R0);	//  Calculate max distance from center
      }				// End loop over dimensions
    R0 *= 1.00001;		// Add some leeway to radius
  }

  //! Build cells of tree adaptively using a top-down approach based on recursion
  void buildCells(Body * bodies, Body * buffer, int begin, int end, Cell * cell, Cells & cells,
		  real_t * X, real_t R, int level = 0, bool direction = false)
  {
    //! Create a tree cell
    cell->BODY = bodies + begin;	// Pointer of first body in cell
    if(direction)
      cell->BODY = buffer + begin;	// Pointer of first body in cell
    cell->NBODY = end - begin;	// Number of bodies in cell
    cell->NCHILD = 0;		// Initialize counter for child cells
    for(int d = 0; d < 3; d++)
      cell->X[d] = X[d];	// Center position of cell
    cell->R = R / (1 << level);	// Cell radius
    //! If cell is a leaf
    if(end - begin <= ncrit)
      {				// If number of bodies is less than threshold
	if(direction)
	  {			//  If direction of data is from bodies to buffer
	    for(int i = begin; i < end; i++)
	      {			//   Loop over bodies in cell
		buffer[i] = bodies[i];	//    Copy bodies source to buffer
	      }			//   End loop over bodies in cell
	  }			//  End if for direction of data
	return;			//  Return without recursion
      }				// End if for number of bodies
    //! Count number of bodies in each octant
    int size[8] = { 0, 0, 0, 0, 0, 0, 0, 0 };
    real_t x[3];		// Coordinates of bodies
    for(int i = begin; i < end; i++)
      {				// Loop over bodies in cell
	for(int d = 0; d < 3; d++)
	  x[d] = bodies[i].X[d];	//  Position of body
	int octant = (x[0] > X[0]) + ((x[1] > X[1]) << 1) + ((x[2] > X[2]) << 2);	// Which octant body belongs to
	size[octant]++;		//  Increment body count in octant
      }				// End loop over bodies in cell
    //! Exclusive scan to get offsets
    int offset = begin;		// Offset of first octant
    int offsets[8], counter[8];	// Offsets and counter for each octant
    for(int i = 0; i < 8; i++)
      {				// Loop over elements
	offsets[i] = offset;	//  Set value
	offset += size[i];	//  Increment offset
	if(size[i])
	  cell->NCHILD++;	//  Increment child cell counter
      }				// End loop over elements
    //! Sort bodies by octant
    for(int i = 0; i < 8; i++)
      counter[i] = offsets[i];	// Copy offsets to counter
    for(int i = begin; i < end; i++)
      {				// Loop over bodies
	for(int d = 0; d < 3; d++)
	  x[d] = bodies[i].X[d];	//  Position of body
	int octant = (x[0] > X[0]) + ((x[1] > X[1]) << 1) + ((x[2] > X[2]) << 2);	// Which octant body belongs to`
	buffer[counter[octant]] = bodies[i];	//  Permute bodies sources out-of-place according to octant
	counter[octant]++;	//  Increment body count in octant
      }				// End loop over bodies
    //! Loop over children and recurse
    real_t Xchild[3];		// Coordinates of children
    cells.resize(cells.size() + cell->NCHILD);	// Resize cell vector
    Cell *child = &cells.back() - cell->NCHILD + 1;	// Pointer for first child cell
    cell->CHILD = child;	// Point to first child cell
    int c = 0;			// Counter for child cells
    for(int i = 0; i < 8; i++)
      {				// Loop over children
	for(int d = 0; d < 3; d++)
	  Xchild[d] = X[d];	//  Initialize center position of child cell
	real_t r = R / (1 << (level + 1));	//  Radius of cells for child's level
	for(int d = 0; d < 3; d++)
	  {			//  Loop over dimensions
	    Xchild[d] += r * (((i & 1 << d) >> d) * 2 - 1);	//   Shift center position to that of child cell
	  }			//  End loop over dimensions
	if(size[i])
	  {			//  If child exists
	    buildCells(buffer, bodies, offsets[i], offsets[i] + size[i],	// Recursive call for each child
		       &child[c], cells, Xchild, R, level + 1, !direction);
	    c++;		//   Increment child cell counter
	  }			//  End if for child
      }				// End loop over children
  }


  void buildCellsSorted(Body * bodies, int begin, int end, Cell * cell, Cells & cells, real_t * X, real_t R, int level = 0)
  {
    //! Create a tree cell
    cell->BODY = bodies + begin;	// Pointer of first body in cell
    cell->NBODY = end - begin;	// Number of bodies in cell
    cell->NCHILD = 0;		// Initialize counter for child cells
    for(int d = 0; d < 3; d++)
      cell->X[d] = X[d];	// Center position of cell
    cell->R = R / (1 << level);	// Cell radius
    //! If cell is a leaf
    if(end - begin <= ncrit)
      return;
    //! Count number of bodies in each octant
    int size[8] = { 0 };
    int offsets[8];		// Offsets and counter for each octant
    real_t x[3];		// Coordinates of bodies

    if(end - begin <= 32)
      {
	for(int i = begin; i < end; i++)
	  {			// Loop over bodies in cell
	    for(int d = 0; d < 3; d++)
	      x[d] = bodies[i].X[d];	//  Position of body
	    int octant = (x[0] >= X[0]) + ((x[1] >= X[1]) << 1) + ((x[2] >= X[2]) << 2);	// Which octant body belongs to
	    size[octant]++;	//  Increment body count in octant
	  }

	//! Exclusive scan to get offsets
	int offset = begin;	// Offset of first octant
	for(int i = 0; i < 8; i++)
	  {			// Loop over elements
	    offsets[i] = offset;	//  Set value
	    offset += size[i];	//  Increment offset
	    if(size[i])
	      cell->NCHILD++;	//  Increment child cell counter
	  }
      }
    else
      {
	offsets[0] = begin;

	for(int i = 1; i < 8; i++)
	  {
	    int left = begin, right = end, mid, octant;
	    while(left < right - 1)
	      {
		mid = (left + right) / 2;
		for(int d = 0; d < 3; d++)
		  x[d] = bodies[mid].X[d];
		octant = (x[0] >= X[0]) + ((x[1] >= X[1]) << 1) + ((x[2] >= X[2]) << 2);
		if(octant >= i)
		  right = mid;
		else
		  left = mid;
	      }
	    offsets[i] = left + 1;
	  }

	for(int i = 0; i < 7; i++)
	  size[i] = offsets[i + 1] - offsets[i];
	size[7] = end - offsets[7];

	for(int i = 0; i < 8; i++)
	  if(size[i])
	    cell->NCHILD++;
      }

    // End loop over bodies
    //! Loop over children and recurse
    real_t Xchild[3];		// Coordinates of children
    cells.resize(cells.size() + cell->NCHILD);	// Resize cell vector
    Cell *child = &cells.back() - cell->NCHILD + 1;	// Pointer for first child cell
    cell->CHILD = child;	// Point to first child cell
    int c = 0;			// Counter for child cells
    for(int i = 0; i < 8; i++)
      {				// Loop over children
	for(int d = 0; d < 3; d++)
	  Xchild[d] = X[d];	//  Initialize center position of child cell
	real_t r = R / (1 << (level + 1));	//  Radius of cells for child's level
	for(int d = 0; d < 3; d++)
	  {			//  Loop over dimensions
	    Xchild[d] += r * (((i & 1 << d) >> d) * 2 - 1);	//   Shift center position to that of child cell
	  }			//  End loop over dimensions
	if(size[i])
	  {			//  If child exists
	    buildCellsSorted(bodies, offsets[i], offsets[i] + size[i],	// Recursive call for each child
			     &child[c], cells, Xchild, R, level + 1);
	    c++;		//   Increment child cell counter
	  }			//  End if for child
      }				// End loop over children
  }


  Cells buildTree(Bodies & bodies)
  {
    real_t R0, X0[3];		// Radius and center root cell
    getBounds(bodies, R0, X0);	// Get bounding box from bodies

    // /*
    unsigned int max_int = ((unsigned int) 1) << (MAX_DEPTH);

#pragma omp parallel for
    for(size_t i = 0; i < bodies.size(); i++)
      {
	real_t x_f = (bodies[i].X[0] - (X0[0] - R0)) / (2.0 * R0);
	real_t y_f = (bodies[i].X[1] - (X0[1] - R0)) / (2.0 * R0);
	real_t z_f = (bodies[i].X[2] - (X0[2] - R0)) / (2.0 * R0);
	bodies[i].Morton[0] = (unsigned int) floor(x_f * (real_t) max_int);
	bodies[i].Morton[1] = (unsigned int) floor(y_f * (real_t) max_int);
	bodies[i].Morton[2] = (unsigned int) floor(z_f * (real_t) max_int);
      }

    //use a openmp sorting to speed things up, with Intel compiler use parallel STL
    //    std::sort(bodies.begin(), bodies.end());

    if(bodies.size() < 1e5)
      std::sort(bodies.begin(), bodies.end());
    else
      omp_par::merge_sort(&bodies[0], &bodies[bodies.size()]);
    //*/
      
    //Bodies buffer = bodies;   // Copy bodies to buffer
    Cells cells(1);		      // Vector of cells
    cells.reserve(bodies.size()/4);	// Reserve memory space
      
    //buildCells(&bodies[0], &buffer[0], 0, bodies.size(), &cells[0], cells, X0, R0);   // Build tree recursively

    buildCellsSorted(&bodies[0], 0, bodies.size(), &cells[0], cells, X0, R0);	// Build tree from sorted bodies
    return cells;		// Return pointer of root cell
  }
}

#endif
