#ifndef traverse_eager_h
#define traverse_eager_h
#include "exafmm.h"
#include "build_tree.h"

namespace exafmm
{
  //! Recursive call to post-order tree traversal for upward pass
  void upwardPass(Cell * Ci)
  {
    for(Cell * Cj = Ci->CHILD; Cj != Ci->CHILD + Ci->NCHILD; Cj++)
      {				// Loop over child cells
#pragma omp task untied if(Cj->NBODY > 100)	//  Start OpenMP task if large enough task
	upwardPass(Cj);		//  Recursive call for child cell
      }				// End loop over child cells
#pragma omp taskwait		// Synchronize OpenMP tasks


    Ci->M.resize(NTERM, 0.0);
    Ci->L.resize(NTERM, 0.0);
    Ci->Pn.resize(P+1,  0.0);
    
    Ci->has_sink = false;

    Ci->R = 1.732 * Ci->R;

    if(Ci->NCHILD == 0)
      P2M(Ci);			// P2M kernel
    else
      M2M(Ci);			// M2M kernel
  }

  //! Upward pass interface
  void upwardPass(Cells & cells)
  {
#pragma omp parallel		// Start OpenMP
#pragma omp single nowait	// Start OpenMP single region with nowait
    upwardPass(&cells[0]);	// Pass root cell to recursive call
  }

  void upwardPass_low(Cell * Ci)
  {
    for(Cell * Cj = Ci->CHILD; Cj != Ci->CHILD + Ci->NCHILD; Cj++)
      {
#pragma omp task untied if(Cj->NBODY > 100)
	upwardPass_low(Cj);
      }
#pragma omp taskwait

    Ci->has_sink = true;

    Ci->M.resize(NTERM, 0.0);
    Ci->L.resize(NTERM, 0.0);
    Ci->Pn.resize(P+1,  0.0);
    
    if(Ci->NCHILD == 0)
      P2M_low(Ci);
    else
      M2M_low(Ci);
  }

  void upwardPass_low(Cells & cells)
  {
#pragma omp parallel
#pragma omp single nowait
    upwardPass_low(&cells[0]);
  }


  void horizontalPassHigh(Cell * Ci, Cell * Cj, bool get_steps)
  {
    for(int d = 0; d < 3; d++)
      dX[d] = Ci->X[d] - Cj->X[d];

    real_t R2 = norm(dX);
    real_t f1 = 10, f2 = 10;

    real_t Ra = Ci->R;
    real_t Rb = Cj->R;
    real_t Ra_p_Rb = (Ra + Rb);

    if(R2 > Ra_p_Rb * Ra_p_Rb)
      {
	real_t R = sqrt(R2);
	real_t Rp = R2;		//std::pow(R, P+2);
	real_t power_Ra = 1;	//std::pow(Ra, P);
	real_t power_Rb = 1;	//std::pow(Rb, P);

	real_t Eba = 0, Eab = 0;

	for(int k = 0; k <= P; k++)
	  {
	    real_t fac = combinator_coef[P][k];
	    Eab += (Ci->Pn[P-k] * power_Rb * fac);
	    Eba += (Cj->Pn[P - k] * power_Ra * fac);
	    power_Ra *= Ra;
	    power_Rb *= Rb;
	    Rp *= R;
	  }

	real_t fac = 8 * fmax(Ra, Rb) / Ra_p_Rb / Rp;
	Eab *= fac;
	f1 = Eab / (force_accuracy * Cj->min_acc);
	Eba *= fac;
	f2 = Eba / (force_accuracy * Ci->min_acc);
      }

    real_t thres = 1;

    //symmetric tree walk, use f2 < thres for unsymmetric walk

    if(f2 < thres && f1 < thres) 
      {
	if((Ci->NBODY <= 2 && Cj->NBODY <= 8) || (Cj->NBODY <= 2 && Ci->NBODY <= 8))
	  {
	    if(get_steps)
	      {
		P2P(Ci, Cj);
	      }
	    else
	      {
		P2P_simple(Ci, Cj);
	      }
	  }
	else	  
	  {
	    if(!get_steps)
	      M2L_rotate(Ci, Cj);	//  M2L kernel
	  }
      }
    else if(Ci->NCHILD == 0 && Cj->NCHILD == 0)
      {				// Else if both cells are leafs
	if(get_steps)
	  {
	    P2P(Ci, Cj);	//  P2P kernel
	  }
	else
	  {
	    P2P_simple(Ci, Cj);
	  }
      }
    else if(Cj->NCHILD == 0 || (Ci->R >= Cj->R && Ci->NCHILD != 0))
      {				// If Cj is leaf or Ci is larger
	for(Cell * ci = Ci->CHILD; ci != Ci->CHILD + Ci->NCHILD; ci++)
	  {			// Loop over Ci's children
#pragma omp task untied if(ci->NBODY > 100)	//   Start OpenMP task if large enough task
	    horizontalPassHigh(ci, Cj, get_steps);	//   Recursive call to target child cells
	  }			//  End loop over Ci's children
      }
    else
      {				// Else if Ci is leaf or Cj is larger
	for(Cell * cj = Cj->CHILD; cj != Cj->CHILD + Cj->NCHILD; cj++)
	  {			// Loop over Cj's children
	    horizontalPassHigh(Ci, cj, get_steps);	//   Recursive call to source child cells
	  }			//  End loop over Cj's children
      }				// End if for leafs and Ci Cj size
  }

  //! Recursive call to dual tree traversal for horizontal pass
  void horizontalPass(Cell * Ci, Cell * Cj, bool get_steps)
  {
    for(int d = 0; d < 3; d++)
      dX[d] = Ci->X[d] - Cj->X[d];	// Distance vector from source to target

    real_t R2 = norm(dX) * theta * theta;	// Scalar distance squared

    if(R2 > (Ci->R + Cj->R) * (Ci->R + Cj->R))
      {				// If distance is far enough
	if((Ci->NBODY <= 2 && Cj->NBODY <= 8) || (Cj->NBODY <= 2 && Ci->NBODY <= 8))
	  {
	    if(get_steps)
	      {
		P2P(Ci, Cj);
	      }
	    else
	      {
		P2P_simple(Ci, Cj);
	      }
	  }
	else
	  {
	    if(!get_steps)
	      M2L_rotate(Ci, Cj);	//  M2L kernel
	  }
      }
    else if(Ci->NCHILD == 0 && Cj->NCHILD == 0)
      {				// Else if both cells are leafs
	if(get_steps)
	  {
	    P2P(Ci, Cj);	//  P2P kernel
	  }
	else
	  {
	    P2P_simple(Ci, Cj);	// simple P2P kernel
	  }
      }
    else if(Cj->NCHILD == 0 || (Ci->R >= Cj->R && Ci->NCHILD != 0))
      {				// If Cj is leaf or Ci is larger
	for(Cell * ci = Ci->CHILD; ci != Ci->CHILD + Ci->NCHILD; ci++)
	  {			// Loop over Ci's children
#pragma omp task untied if(ci->NBODY > 100)	//   Start OpenMP task if large enough task
	    horizontalPass(ci, Cj, get_steps);	//   Recursive call to target child cells
	  }			//  End loop over Ci's children
      }
    else
      {				// Else if Ci is leaf or Cj is larger
	for(Cell * cj = Cj->CHILD; cj != Cj->CHILD + Cj->NCHILD; cj++)
	  {			// Loop over Cj's children
	    horizontalPass(Ci, cj, get_steps);	//   Recursive call to source child cells
	  }			//  End loop over Cj's children
      }				// End if for leafs and Ci Cj size
  }

  //! Horizontal pass interface
  void horizontalPass(Cells & icells, Cells & jcells, bool high_force, bool get_steps)
  {
#pragma omp parallel		// Start OpenMP
#pragma omp single nowait
    {				// Start OpenMP single region with nowait
      if(!high_force)
	horizontalPass(&icells[0], &jcells[0], get_steps);	// Pass root cell to recursive call
      else
	horizontalPassHigh(&icells[0], &jcells[0], get_steps);	// Pass root cell to recursive call
    }
  }


  //! Recursive call to dual tree traversal for horizontal pass
  void horizontalPass_low(Cell * Ci, Cell * Cj)
  {
    for(int d = 0; d < 3; d++)
      dX[d] = Ci->X[d] - Cj->X[d];
    real_t R2 = norm(dX);
    //use theta=1; notice C->R is the side length of a cell, therefore R should multiply by sqrt(3)
    if(R2 > (Ci->R + Cj->R) * (Ci->R + Cj->R) * 3)
      {
	M2L_low(Ci, Cj);
      }
    else if(Ci->NCHILD == 0 && Cj->NCHILD == 0)
      {
	P2P_low(Ci, Cj);
      }
    else if(Cj->NCHILD == 0 || (Ci->R >= Cj->R && Ci->NCHILD != 0))
      {
	for(Cell * ci = Ci->CHILD; ci != Ci->CHILD + Ci->NCHILD; ci++)
	  {
#pragma omp task untied if(ci->NBODY > 1000)
	    horizontalPass_low(ci, Cj);
	  }
      }
    else
      {
	for(Cell * cj = Cj->CHILD; cj != Cj->CHILD + Cj->NCHILD; cj++)
	  {
	    horizontalPass_low(Ci, cj);
	  }
      }
  }

  //! Horizontal pass interface
  void horizontalPass_low(Cells & icells, Cells & jcells)
  {
#pragma omp parallel
#pragma omp single nowait
    {
      horizontalPass_low(&icells[0], &jcells[0]);
    }
  }


  void directPass(Cell * Ci, Cell * Cj, bool get_steps)
  {
    if(Ci->NCHILD == 0 && Cj->NCHILD == 0)
      {				// Else if both cells are leafs
	if(get_steps)
	  {
	    P2P(Ci, Cj);
	  }
	else
	  {
	    P2P_simple(Ci, Cj);
	  }			//  P2P kernel
      }
    else if(Cj->NCHILD == 0 || (Ci->R >= Cj->R && Ci->NCHILD != 0))
      {				// If Cj is leaf or Ci is larger
	for(Cell * ci = Ci->CHILD; ci != Ci->CHILD + Ci->NCHILD; ci++)
	  {			// Loop over Ci's children
#pragma omp task untied if(ci->NBODY > 100)	//   Start OpenMP task if large enough task
	    directPass(ci, Cj, get_steps);	//   Recursive call to target child cells
	  }			//  End loop over Ci's children
      }
    else
      {				// Else if Ci is leaf or Cj is larger
	for(Cell * cj = Cj->CHILD; cj != Cj->CHILD + Cj->NCHILD; cj++)
	  {			// Loop over Cj's children
	    directPass(Ci, cj, get_steps);	//   Recursive call to source child cells
	  }			//  End loop over Cj's children
      }				// End if for leafs and Ci Cj size
  }

  //! direct pass interface
  void directPass(Cells & icells, Cells & jcells, bool get_steps)
  {
#pragma omp parallel		// Start OpenMP
#pragma omp single nowait
    {				// Start OpenMP single region with nowait
      directPass(&icells[0], &jcells[0], get_steps);	// Pass root cell to recursive call
    }
  }

  //! Recursive call to pre-order tree traversal for downward pass
  void downwardPass_low(Cell * Cj)
  {
    if(Cj->NCHILD == 0)
      L2P_low(Cj);		// L2P kernel
    else
      L2L_low(Cj);
    for(Cell * Ci = Cj->CHILD; Ci != Cj->CHILD + Cj->NCHILD; Ci++)
      {
#pragma omp task untied if(Ci->NBODY > 100)
	downwardPass_low(Ci);
      }
#pragma omp taskwait
  }

  //! Downward pass interface
  void downwardPass_low(Cells & cells)
  {
#pragma omp parallel
#pragma omp single nowait
    downwardPass_low(&cells[0]);
  }

  //! Recursive call to pre-order tree traversal for downward pass
  void downwardPass(Cell * Cj)
  {
    if(Cj->NCHILD == 0)
      L2P(Cj);			// L2P kernel
    else
      L2L(Cj);
    for(Cell * Ci = Cj->CHILD; Ci != Cj->CHILD + Cj->NCHILD; Ci++)
      {				// Loop over child cells
#pragma omp task untied if(Ci->NBODY > 100)	//  Start OpenMP task if large enough task
	downwardPass(Ci);	//  Recursive call for child cell
      }				// End loop over chlid cells
#pragma omp taskwait		// Synchronize OpenMP tasks
  }

  //! Downward pass interface
  void downwardPass(Cells & cells)
  {
#pragma omp parallel		// Start OpenMP
#pragma omp single nowait	// Start OpenMP single region with nowait
    downwardPass(&cells[0]);	// Pass root cell to recursive call
  }

  //! Direct summation
  void direct(Bodies & bodies, Bodies & jbodies, bool get_steps)
  {
    Cells cells = buildTree(bodies);
    Cells jcells = buildTree(jbodies);
//    for(size_t b = 0; b < cells.size(); b++)
//      {
//	cells[b].has_sink = true;
//      }
    upwardPass(cells);
    upwardPass(jcells);
    directPass(cells, jcells, get_steps);

    //Cells cells(2);                                             // Define a pair of cells to pass to P2P kernel
    //Cell * Ci = &cells[0];                                      // Allocate single target
    //Cell * Cj = &cells[1];                                      // Allocate single source
    //Ci->BODY = &bodies[0];                                      // Iterator of first target body
    //Ci->NBODY = bodies.size();                                  // Number of target bodies
    //Cj->BODY = &jbodies[0];                                     // Iterator of first source body
    //Cj->NBODY = jbodies.size();                                 // Number of source bodies
    //Ci->has_sink = true;
    //P2P_OMP(Ci, Cj);                                                // Evaluate P2P kenrel
  }
}
#endif
