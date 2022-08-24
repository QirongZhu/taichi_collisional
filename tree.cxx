#include "tree.h"

int main()
{
  //allocate bodies
  Particles particles(100);
  
  //allocate cells  
  Cells cells(particles.size() - 1);
  
  return 0;
}
