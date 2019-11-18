#include "exafmm.h"

#include "vectorclass/vectorclass.h"

namespace exafmm
{

  void swap_x_z_two_arrays(real_t * arr1, real_t * arr2)
  {
    Vec4d arrold[NTERM];
    Vec4d arr[NTERM];

    for(int i = 0; i < NTERM; i++)
      {
	Vec4d temp(arr1[i], arr2[i], 0.0, 0.0);
	  arrold[i] = temp;
      }

    arr[0] = arrold[0];
    arr[1] = arrold[1];
    arr[2] = arrold[3];
    arr[3] = arrold[2];

    real_t coeff;

#if EXPANSION > 1
    coeff = 1 / 2.0;
    arr[4] = (4 * arrold[5]) * coeff;
    arr[5] = (arrold[4]) * coeff;
    arr[6] = (-arrold[6] + arrold[8]) * coeff;
    arr[7] = (2 * arrold[7]) * coeff;
    arr[8] = (3 * arrold[6] + arrold[8]) * coeff;
#endif

#include "multipole_rotation_matrix/rotate_order_3.txt"
#include "multipole_rotation_matrix/rotate_order_4.txt"
#include "multipole_rotation_matrix/rotate_order_5.txt"
#include "multipole_rotation_matrix/rotate_order_6.txt"
#include "multipole_rotation_matrix/rotate_order_7.txt"
#include "multipole_rotation_matrix/rotate_order_8.txt"
#include "multipole_rotation_matrix/rotate_order_9.txt"
#include "multipole_rotation_matrix/rotate_order_10.txt"
#include "multipole_rotation_matrix/rotate_order_11.txt"
#include "multipole_rotation_matrix/rotate_order_12.txt"
#include "multipole_rotation_matrix/rotate_order_13.txt"
#include "multipole_rotation_matrix/rotate_order_14.txt"
#include "multipole_rotation_matrix/rotate_order_15.txt"
#include "multipole_rotation_matrix/rotate_order_16.txt"
#include "multipole_rotation_matrix/rotate_order_17.txt"
#include "multipole_rotation_matrix/rotate_order_18.txt"
#include "multipole_rotation_matrix/rotate_order_19.txt"
#include "multipole_rotation_matrix/rotate_order_20.txt"
#include "multipole_rotation_matrix/rotate_order_21.txt"
#include "multipole_rotation_matrix/rotate_order_22.txt"
#include "multipole_rotation_matrix/rotate_order_23.txt"
#include "multipole_rotation_matrix/rotate_order_24.txt"
#include "multipole_rotation_matrix/rotate_order_25.txt"
#include "multipole_rotation_matrix/rotate_order_26.txt"
#include "multipole_rotation_matrix/rotate_order_27.txt"
#include "multipole_rotation_matrix/rotate_order_28.txt"
#include "multipole_rotation_matrix/rotate_order_29.txt"
#include "multipole_rotation_matrix/rotate_order_30.txt"

    for(int i = 0; i < NTERM; i++)
      {
	Vec4d temp = arr[i];
	arr1[i] = temp[0];
	arr2[i] = temp[1];
      }

    
  }
}
