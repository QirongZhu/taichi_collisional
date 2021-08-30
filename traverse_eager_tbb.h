#ifndef traverse_eager_h
#define traverse_eager_h
#include "exafmm.h"
#include "build_tree.h"

namespace exafmm
{
  //! Recursive call to post-order tree traversal for upward pass
  void upwardPass(Cell * Ci)
  {
    tbb::task_group tg;
	  
    for(Cell * Cj = Ci->CHILD; Cj != Ci->CHILD + Ci->NCHILD; Cj++) {
      if(Cj->NBODY > ncrit) {
	tg.run([=]{ upwardPass(Cj);});
      }
      else{
	upwardPass(Cj);
      }
    }

    tg.wait();
    
    Ci->has_sink = false;
    Ci->min_acc  = HUGE;

    Ci->R = 1.732 * Ci->R;

    if(Ci->NCHILD == 0){
      P2M(Ci);			// P2M kernel
    }
    else {
      M2M(Ci);			// M2M kernel
    }
  }

  //! Upward pass interface
  void upwardPass(Cells & cells)
  {
    upwardPass(&cells[0]);	// Pass root cell to recursive call
  }

  void upwardPass_low(Cell * Ci)
  {
    tbb::task_group tg;
      
    for(Cell * Cj = Ci->CHILD; Cj != Ci->CHILD + Ci->NCHILD; Cj++) {
      if(Cj->NBODY > ncrit) {
	tg.run([=] {upwardPass_low(Cj);});
      }else{
	upwardPass_low(Cj);
      }
    }
      
    tg.wait();

    Ci->has_sink = true;

    if(Ci->NCHILD == 0)
      P2M_low(Ci);
    else
      M2M_low(Ci);
  }

  void upwardPass_low(Cells & cells)
  {
    upwardPass_low(&cells[0]);
  }

  void horizontalPassHigh(Cell * Ci, Cell * Cj)
  {
    if(!Ci->has_sink || std::abs(Cj->M[0]) < 1e-16) {
      return;
    }
      
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
	real_t Rp = R2;		    //std::pow(R, P+2);
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
	M2L_rotate(Ci, Cj);	//  M2L kernel
      }
    else if(Ci->NCHILD == 0 && Cj->NCHILD == 0)
      {				// Else if both cells are leafs
	P2P(Ci, Cj);	//  P2P kernel
      }
    else if(Cj->NCHILD == 0 || (Ci->R >= Cj->R && Ci->NCHILD != 0))
      {				// If Cj is leaf or Ci is larger

	tbb::task_group tg;
	
	for(Cell * ci = Ci->CHILD; ci != Ci->CHILD + Ci->NCHILD; ci++)
	  {			// Loop over Ci's children
	    if(ci->NBODY > USETBB){
	      tg.run([=]{horizontalPassHigh(ci, Cj);});
	    }
	    else{
	      horizontalPassHigh(ci, Cj);
	    }			    
	  }			//  End loop over Ci's children	

	tg.wait();
	
      }
    else
      {				// Else if Ci is leaf or Cj is larger
	for(Cell * cj = Cj->CHILD; cj != Cj->CHILD + Cj->NCHILD; cj++)
	  {			// Loop over Cj's children
	    horizontalPassHigh(Ci, cj);	//Recursive call to source child cells
	  }			//  End loop over Cj's children
      }				// End if for leafs and Ci Cj size
  }

  //! Recursive call to dual tree traversal for horizontal pass
  void horizontalPass(Cell * Ci, Cell * Cj)
  {
    for(int d = 0; d < 3; d++)
      dX[d] = Ci->X[d] - Cj->X[d];  // Distance vector from source to target

    real_t R2 = norm(dX) * theta * theta;	// Scalar distance squared

    if(R2 > (Ci->R + Cj->R) * (Ci->R + Cj->R))
      {				// If distance is far enough
	M2L_rotate(Ci, Cj);	//  M2L kernel
      }
    else if(Ci->NCHILD == 0 && Cj->NCHILD == 0)
      {				// Else if both cells are leafs
	P2P(Ci, Cj);	//  P2P kernel
      }
    else if(Cj->NCHILD == 0 || (Ci->R >= Cj->R && Ci->NCHILD != 0))
      {				// If Cj is leaf or Ci is larger
        
	tbb::task_group tg;

	for(Cell * ci = Ci->CHILD; ci != Ci->CHILD + Ci->NCHILD; ci++)
	  {			// Loop over Ci's children
	    if(ci->NBODY > USETBB){
              tg.run([=] {horizontalPass(ci, Cj);});
	    }else{
              horizontalPass(ci, Cj);
	    }
	  }
          
	tg.wait();
          
      }
    else
      {				// Else if Ci is leaf or Cj is larger
	for(Cell * cj = Cj->CHILD; cj != Cj->CHILD + Cj->NCHILD; cj++)
	  {			// Loop over Cj's children
	    horizontalPass(Ci, cj); // Recursive call to source child cells
	  }			//  End loop over Cj's children
      }				// End if for leafs and Ci Cj size
  }

  //! Horizontal pass interface
  void horizontalPass(Cells & icells, Cells & jcells, bool high_force)
  {
    if(!high_force)
      horizontalPass(&icells[0], &jcells[0]);//Pass root cell to recursive call
    else
      horizontalPassHigh(&icells[0], &jcells[0]); 
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
    else
      {
	if(Cj->NCHILD == 0 || (Ci->R >= Cj->R && Ci->NCHILD != 0))
	  {
	    tbb::task_group tg;

	    for(Cell * ci = Ci->CHILD; ci != Ci->CHILD + Ci->NCHILD; ci++)
	      {
		if(ci->NBODY > USETBB) {
          tg.run([=] {horizontalPass_low(ci, Cj);});
		}
		else
		  horizontalPass_low(ci, Cj);
	      }
	    tg.wait();
	  }
	else
	  {
	    for(Cell * cj = Cj->CHILD; cj != Cj->CHILD + Cj->NCHILD; cj++)
	      {
		horizontalPass_low(Ci, cj);
	      }
	  }
      }
  }

  //! Horizontal pass interface
  void horizontalPass_low(Cells & icells, Cells & jcells)
  {
    horizontalPass_low(&icells[0], &jcells[0]);
  }


  void directPass(Cell * Ci, Cell * Cj)
  {
    if(Ci->NCHILD == 0 && Cj->NCHILD == 0) {// Else if both cells are leafs
      P2P(Ci, Cj);
    }
    else if(Cj->NCHILD == 0 || (Ci->R >= Cj->R && Ci->NCHILD != 0)) {// If Cj is leaf or Ci is larger

      tbb::task_group tg;

      for(Cell * ci = Ci->CHILD; ci != Ci->CHILD + Ci->NCHILD; ci++) {// Loop over Ci's children
	if(ci->NBODY > ncrit/4) {
	  tg.run([=]{directPass(ci, Cj);});
	}
	else{
	  directPass(ci, Cj);
	}
      }
          
      tg.wait();
          
    }
    else {// Else if Ci is leaf or Cj is larger
      for(Cell * cj = Cj->CHILD; cj != Cj->CHILD + Cj->NCHILD; cj++) {// Loop over Cj's children
	directPass(Ci, cj);	//   Recursive call to source child cells
      }			//  End loop over Cj's children
    }				// End if for leafs and Ci Cj size
  }

  //! direct pass interface
  void directPass(Cells & icells, Cells & jcells) {
    directPass(&icells[0], &jcells[0]);	// Pass root cell to recursive call
  }

  //! Recursive call to pre-order tree traversal for downward pass
  void downwardPass_low(Cell * Cj)
  {
    if(Cj->NCHILD == 0)
      L2P_low(Cj);		// L2P kernel
    else
      L2L_low(Cj);

    tbb::task_group tg;
        
    for(Cell * Ci = Cj->CHILD; Ci != Cj->CHILD + Cj->NCHILD; Ci++) {
      if(Ci->NBODY > ncrit) {
	tg.run([=]{ downwardPass_low(Ci);});
      }else{
	downwardPass_low(Ci);
      }
    }

    tg.wait();
        
  }

  //! Downward pass interface
  void downwardPass_low(Cells & cells) {
    downwardPass_low(&cells[0]);
  }

  //! Recursive call to pre-order tree traversal for downward pass
  void downwardPass(Cell * Cj)
  {
    if(Cj->NCHILD == 0)
      L2P(Cj);			// L2P kernel
    else
      L2L(Cj);

    tbb::task_group tg;
    
    for(Cell * Ci = Cj->CHILD; Ci != Cj->CHILD + Cj->NCHILD; Ci++) {// Loop over child cells
      if(Ci->NBODY > ncrit) {
	tg.run([=]{ downwardPass(Ci);});
      }else{
	downwardPass(Ci);
      }
    }				// End loop over chlid cells
    tg.wait();
  }

  //! Downward pass interface
  void downwardPass(Cells & cells)
  {
    downwardPass(&cells[0]);	// Pass root cell to recursive call
  }

  //! Direct summation
  void direct(Bodies & bodies)
  {
    Cells cells = buildTree(bodies);
            
    for (size_t i=0; i<cells.size(); i++) {
      cells[i].has_sink = true;
    }
    
    directPass(cells, cells);// Evaluate P2P kenrel
  }
}
#endif
