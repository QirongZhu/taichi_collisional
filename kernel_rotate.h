#include "exafmm.h"

namespace exafmm
{

  void swap_x_z(real_t * arr)
  {
    real_t arrold[NTERM];
      memcpy(arrold, arr, NTERM * sizeof(real_t));

      arr[2] = arrold[3];
      arr[3] = arrold[2];

    real_t coeff;

#if EXPANSION > 1
      coeff = 1 / 2.0;
      arr[4] = (4 * arrold[5]) * coeff;
      arr[5] = (arrold[4]) * coeff;
      arr[6] = (-arrold[6] + arrold[8]) * coeff;
      /*arr[7] = (2*arrold[7])*coeff;*/
      arr[8] = (3 * arrold[6] + arrold[8]) * coeff;
#endif

#if EXPANSION > 2
      coeff = 1 / 4.0;
      arr[9] = (arrold[9] + 15 * arrold[11]) * coeff;
      /*arr[10] = (4*arrold[10])*coeff; */
      arr[11] = (arrold[9] - arrold[11]) * coeff;
      arr[12] = (-3 * arrold[13] + arrold[15]) * coeff;
      arr[13] = (-2 * arrold[12] + 2 * arrold[14]) * coeff;
      arr[14] = (5 * arrold[13] + arrold[15]) * coeff;
      arr[15] = (10 * arrold[12] + 6 * arrold[14]) * coeff;
#endif
#if EXPANSION > 3
      coeff = 1 / 8.0;
      arr[16] = (8 * arrold[17] + 56 * arrold[19]) * coeff;
      arr[17] = (arrold[16] + 14 * arrold[18]) * coeff;
      arr[18] = (4 * arrold[17] - 4 * arrold[19]) * coeff;
      arr[19] = (arrold[16] - 2 * arrold[18]) * coeff;
      arr[20] = (3 * arrold[20] - 4 * arrold[22] + arrold[24]) * coeff;
      arr[21] = (-6 * arrold[21] + 2 * arrold[23]) * coeff;
      arr[22] = (-5 * arrold[20] + 4 * arrold[22] + arrold[24]) * coeff;
      arr[23] = (14 * arrold[21] + 6 * arrold[23]) * coeff;
      arr[24] = (35 * arrold[20] + 28 * arrold[22] + arrold[24]) * coeff;
#endif
#if EXPANSION > 4
      coeff = 1 / 16.0;
      arr[25] = (arrold[25] + 45 * arrold[27] + 210 * arrold[29]) * coeff;
      arr[26] = (8 * arrold[26] + 48 * arrold[28]) * coeff;
      arr[27] = (arrold[25] + 13 * arrold[27] - 14 * arrold[29]) * coeff;
      arr[28] = (4 * arrold[26] - 8 * arrold[28]) * coeff;
      arr[29] = (arrold[25] - 3 * arrold[27] + 2 * arrold[29]) * coeff;
      arr[30] = (10 * arrold[31] - 5 * arrold[33] + arrold[35]) * coeff;
      arr[31] = (6 * arrold[30] - 8 * arrold[32] + 2 * arrold[34]) * coeff;
      arr[32] = (-14 * arrold[31] + 3 * arrold[33] + arrold[35]) * coeff;
      arr[33] = (-14 * arrold[30] + 8 * arrold[32] + 6 * arrold[34]) * coeff;
      arr[34] = (42 * arrold[31] + 27 * arrold[33] + arrold[35]) * coeff;
      arr[35] = (126 * arrold[30] + 120 * arrold[32] + 10 * arrold[34]) * coeff;
#endif
#if EXPANSION > 5
      coeff = 1 / 32.0;
      arr[36] = (12 * arrold[37] + 220 * arrold[39] + 792 * arrold[41]) * coeff;
      arr[37] = (arrold[36] + 44 * arrold[38] + 165 * arrold[40]) * coeff;
      arr[38] = (8 * arrold[37] + 40 * arrold[39] - 48 * arrold[41]) * coeff;
      arr[39] = (arrold[36] + 12 * arrold[38] - 27 * arrold[40]) * coeff;
      arr[40] = (4 * arrold[37] - 12 * arrold[39] + 8 * arrold[41]) * coeff;
      arr[41] = (arrold[36] - 4 * arrold[38] + 5 * arrold[40]) * coeff;
      arr[42] = (-10 * arrold[42] + 15 * arrold[44] - 6 * arrold[46] + arrold[48]) * coeff;
      arr[43] = (20 * arrold[43] - 10 * arrold[45] + 2 * arrold[47]) * coeff;
      arr[44] = (14 * arrold[42] - 17 * arrold[44] + 2 * arrold[46] + arrold[48]) * coeff;
      arr[45] = (-36 * arrold[43] + 2 * arrold[45] + 6 * arrold[47]) * coeff;
      arr[46] = (-42 * arrold[42] + 15 * arrold[44] + 26 * arrold[46] + arrold[48]) * coeff;
      arr[47] = (132 * arrold[43] + 110 * arrold[45] + 10 * arrold[47]) * coeff;
      arr[48] = (462 * arrold[42] + 495 * arrold[44] + 66 * arrold[46] + arrold[48]) * coeff;
#endif
#if EXPANSION > 6
      coeff = 1 / 64.0;
      arr[49] = (arrold[49] + 91 * arrold[51] + 1001 * arrold[53] + 3003 * arrold[55]) * coeff;
      arr[50] = (12 * arrold[50] + 208 * arrold[52] + 572 * arrold[54]) * coeff;
      arr[51] = (arrold[49] + 43 * arrold[51] + 121 * arrold[53] - 165 * arrold[55]) * coeff;
      arr[52] = (8 * arrold[50] + 32 * arrold[52] - 88 * arrold[54]) * coeff;
      arr[53] = (arrold[49] + 11 * arrold[51] - 39 * arrold[53] + 27 * arrold[55]) * coeff;
      arr[54] = (4 * arrold[50] - 16 * arrold[52] + 20 * arrold[54]) * coeff;
      arr[55] = (arrold[49] - 5 * arrold[51] + 9 * arrold[53] - 5 * arrold[55]) * coeff;
      arr[56] = (-35 * arrold[57] + 21 * arrold[59] - 7 * arrold[61] + arrold[63]) * coeff;
      arr[57] = (-20 * arrold[56] + 30 * arrold[58] - 12 * arrold[60] + 2 * arrold[62]) * coeff;
      arr[58] = (45 * arrold[57] - 19 * arrold[59] + arrold[61] + arrold[63]) * coeff;
      arr[59] = (36 * arrold[56] - 38 * arrold[58] - 4 * arrold[60] + 6 * arrold[62]) * coeff;
      arr[60] = (-99 * arrold[57] - 11 * arrold[59] + 25 * arrold[61] + arrold[63]) * coeff;
      arr[61] = (-132 * arrold[56] + 22 * arrold[58] + 100 * arrold[60] + 10 * arrold[62]) * coeff;
      arr[62] = (429 * arrold[57] + 429 * arrold[59] + 65 * arrold[61] + arrold[63]) * coeff;
      arr[63] = (1716 * arrold[56] + 2002 * arrold[58] + 364 * arrold[60] + 14 * arrold[62]) * coeff;
#endif
#if EXPANSION > 7
      coeff = 1 / 128.0;
      arr[64] = (16 * arrold[65] + 560 * arrold[67] + 4368 * arrold[69] + 11440 * arrold[71]) * coeff;
      arr[65] = (arrold[64] + 90 * arrold[66] + 910 * arrold[68] + 2002 * arrold[70]) * coeff;
      arr[66] = (12 * arrold[65] + 196 * arrold[67] + 364 * arrold[69] - 572 * arrold[71]) * coeff;
      arr[67] = (arrold[64] + 42 * arrold[66] + 78 * arrold[68] - 286 * arrold[70]) * coeff;
      arr[68] = (8 * arrold[65] + 24 * arrold[67] - 120 * arrold[69] + 88 * arrold[71]) * coeff;
      arr[69] = (arrold[64] + 10 * arrold[66] - 50 * arrold[68] + 66 * arrold[70]) * coeff;
      arr[70] = (4 * arrold[65] - 20 * arrold[67] + 36 * arrold[69] - 20 * arrold[71]) * coeff;
      arr[71] = (arrold[64] - 6 * arrold[66] + 14 * arrold[68] - 14 * arrold[70]) * coeff;
      arr[72] = (35 * arrold[72] - 56 * arrold[74] + 28 * arrold[76] - 8 * arrold[78] + arrold[80]) * coeff;
      arr[73] = (-70 * arrold[73] + 42 * arrold[75] - 14 * arrold[77] + 2 * arrold[79]) * coeff;
      arr[74] = (-45 * arrold[72] + 64 * arrold[74] - 20 * arrold[76] + arrold[80]) * coeff;
      arr[75] = (110 * arrold[73] - 34 * arrold[75] - 10 * arrold[77] + 6 * arrold[79]) * coeff;
      arr[76] = (99 * arrold[72] - 88 * arrold[74] - 36 * arrold[76] + 24 * arrold[78] + arrold[80]) * coeff;
      arr[77] = (-286 * arrold[73] - 78 * arrold[75] + 90 * arrold[77] + 10 * arrold[79]) * coeff;
      arr[78] = (-429 * arrold[72] + 364 * arrold[76] + 64 * arrold[78] + arrold[80]) * coeff;
      arr[79] = (1430 * arrold[73] + 1638 * arrold[75] + 350 * arrold[77] + 14 * arrold[79]) * coeff;
      arr[80] =
      (6435 * arrold[72] + 8008 * arrold[74] + 1820 * arrold[76] + 120 * arrold[78] + arrold[80]) * coeff;
#endif
#if EXPANSION > 8
      coeff = 1 / 256.0;
      arr[81] =
      (arrold[81] + 153 * arrold[83] + 3060 * arrold[85] + 18564 * arrold[87] + 43758 * arrold[89]) * coeff;
      arr[82] = (16 * arrold[82] + 544 * arrold[84] + 3808 * arrold[86] + 7072 * arrold[88]) * coeff;
      arr[83] =
      (arrold[81] + 89 * arrold[83] + 820 * arrold[85] + 1092 * arrold[87] - 2002 * arrold[89]) * coeff;
      arr[84] = (12 * arrold[82] + 184 * arrold[84] + 168 * arrold[86] - 936 * arrold[88]) * coeff;
      arr[85] =
      (arrold[81] + 41 * arrold[83] + 36 * arrold[85] - 364 * arrold[87] + 286 * arrold[89]) * coeff;
      arr[86] = (8 * arrold[82] + 16 * arrold[84] - 144 * arrold[86] + 208 * arrold[88]) * coeff;
      arr[87] = (arrold[81] + 9 * arrold[83] - 60 * arrold[85] + 116 * arrold[87] - 66 * arrold[89]) * coeff;
      arr[88] = (4 * arrold[82] - 24 * arrold[84] + 56 * arrold[86] - 56 * arrold[88]) * coeff;
      arr[89] = (arrold[81] - 7 * arrold[83] + 20 * arrold[85] - 28 * arrold[87] + 14 * arrold[89]) * coeff;
      arr[90] = (126 * arrold[91] - 84 * arrold[93] + 36 * arrold[95] - 9 * arrold[97] + arrold[99]) * coeff;
      arr[91] =
      (70 * arrold[90] - 112 * arrold[92] + 56 * arrold[94] - 16 * arrold[96] + 2 * arrold[98]) * coeff;
      arr[92] = (-154 * arrold[91] + 84 * arrold[93] - 20 * arrold[95] - arrold[97] + arrold[99]) * coeff;
      arr[93] =
      (-110 * arrold[90] + 144 * arrold[92] - 24 * arrold[94] - 16 * arrold[96] + 6 * arrold[98]) * coeff;
      arr[94] = (286 * arrold[91] - 52 * arrold[93] - 60 * arrold[95] + 23 * arrold[97] + arrold[99]) * coeff;
      arr[95] =
      (286 * arrold[90] - 208 * arrold[92] - 168 * arrold[94] + 80 * arrold[96] + 10 * arrold[98]) * coeff;
      arr[96] =
      (-858 * arrold[91] - 364 * arrold[93] + 300 * arrold[95] + 63 * arrold[97] + arrold[99]) * coeff;
      arr[97] =
      (-1430 * arrold[90] - 208 * arrold[92] + 1288 * arrold[94] + 336 * arrold[96] +
       14 * arrold[98]) * coeff;
      arr[98] =
      (4862 * arrold[91] + 6188 * arrold[93] + 1700 * arrold[95] + 119 * arrold[97] + arrold[99]) * coeff;
      arr[99] =
      (24310 * arrold[90] + 31824 * arrold[92] + 8568 * arrold[94] + 816 * arrold[96] +
       18 * arrold[98]) * coeff;
#endif
#if EXPANSION > 9
      coeff = 1 / 512.0;
      arr[100] =
      (20 * arrold[101] + 1140 * arrold[103] + 15504 * arrold[105] + 77520 * arrold[107] +
       167960 * arrold[109]) * coeff;
      arr[101] =
      (arrold[100] + 152 * arrold[102] + 2907 * arrold[104] + 15504 * arrold[106] +
       25194 * arrold[108]) * coeff;
      arr[102] =
      (16 * arrold[101] + 528 * arrold[103] + 3264 * arrold[105] + 3264 * arrold[107] -
       7072 * arrold[109]) * coeff;
      arr[103] =
      (arrold[100] + 88 * arrold[102] + 731 * arrold[104] + 272 * arrold[106] - 3094 * arrold[108]) * coeff;
      arr[104] =
      (12 * arrold[101] + 172 * arrold[103] - 16 * arrold[105] - 1104 * arrold[107] +
       936 * arrold[109]) * coeff;
      arr[105] =
      (arrold[100] + 40 * arrold[102] - 5 * arrold[104] - 400 * arrold[106] + 650 * arrold[108]) * coeff;
      arr[106] =
      (8 * arrold[101] + 8 * arrold[103] - 160 * arrold[105] + 352 * arrold[107] - 208 * arrold[109]) * coeff;
      arr[107] =
      (arrold[100] + 8 * arrold[102] - 69 * arrold[104] + 176 * arrold[106] - 182 * arrold[108]) * coeff;
      arr[108] =
      (4 * arrold[101] - 28 * arrold[103] + 80 * arrold[105] - 112 * arrold[107] + 56 * arrold[109]) * coeff;
      arr[109] =
      (arrold[100] - 8 * arrold[102] + 27 * arrold[104] - 48 * arrold[106] + 42 * arrold[108]) * coeff;
      arr[110] =
      (-126 * arrold[110] + 210 * arrold[112] - 120 * arrold[114] + 45 * arrold[116] - 10 * arrold[118] +
       arrold[120]) * coeff;
      arr[111] =
      (252 * arrold[111] - 168 * arrold[113] + 72 * arrold[115] - 18 * arrold[117] + 2 * arrold[119]) * coeff;
      arr[112] =
      (154 * arrold[110] - 238 * arrold[112] + 104 * arrold[114] - 19 * arrold[116] - 2 * arrold[118] +
       arrold[120]) * coeff;
      arr[113] =
      (-364 * arrold[111] + 168 * arrold[113] - 8 * arrold[115] - 22 * arrold[117] + 6 * arrold[119]) * coeff;
      arr[114] =
      (-286 * arrold[110] + 338 * arrold[112] + 8 * arrold[114] - 83 * arrold[116] + 22 * arrold[118] +
       arrold[120]) * coeff;
      arr[115] =
      (780 * arrold[111] - 40 * arrold[113] - 248 * arrold[115] + 70 * arrold[117] +
       10 * arrold[119]) * coeff;
      arr[116] =
      (858 * arrold[110] - 494 * arrold[112] - 664 * arrold[114] + 237 * arrold[116] + 62 * arrold[118] +
       arrold[120]) * coeff;
      arr[117] =
      (-2652 * arrold[111] - 1496 * arrold[113] + 952 * arrold[115] + 322 * arrold[117] +
       14 * arrold[119]) * coeff;
      arr[118] =
      (-4862 * arrold[110] - 1326 * arrold[112] + 4488 * arrold[114] + 1581 * arrold[116] +
       118 * arrold[118] + arrold[120]) * coeff;
      arr[119] =
      (16796 * arrold[111] + 23256 * arrold[113] + 7752 * arrold[115] + 798 * arrold[117] +
       18 * arrold[119]) * coeff;
      arr[120] =
      (92378 * arrold[110] + 125970 * arrold[112] + 38760 * arrold[114] + 4845 * arrold[116] +
       190 * arrold[118] + arrold[120]) * coeff;
#endif

#if EXPANSION > 10
      coeff = 1 / 1024.0;
      arr[121] =
      (arrold[121] + 231 * arrold[123] + 7315 * arrold[125] + 74613 * arrold[127] + 319770 * arrold[129] +
       646646 * arrold[131]) * coeff;
      arr[122] =
      (20 * arrold[122] + 1120 * arrold[124] + 14364 * arrold[126] + 62016 * arrold[128] +
       90440 * arrold[130]) * coeff;
      arr[123] =
      (arrold[121] + 151 * arrold[123] + 2755 * arrold[125] + 12597 * arrold[127] + 9690 * arrold[129] -
       25194 * arrold[131]) * coeff;
      arr[124] = (16 * arrold[122] + 512 * arrold[124] + 2736 * arrold[126] - 10336 * arrold[130]) * coeff;
      arr[125] =
      (arrold[121] + 87 * arrold[123] + 643 * arrold[125] - 459 * arrold[127] - 3366 * arrold[129] +
       3094 * arrold[131]) * coeff;
      arr[126] =
      (12 * arrold[122] + 160 * arrold[124] - 188 * arrold[126] - 1088 * arrold[128] +
       2040 * arrold[130]) * coeff;
      arr[127] =
      (arrold[121] + 39 * arrold[123] - 45 * arrold[125] - 395 * arrold[127] + 1050 * arrold[129] -
       650 * arrold[131]) * coeff;
      arr[128] = (8 * arrold[122] - 168 * arrold[126] + 512 * arrold[128] - 560 * arrold[130]) * coeff;
      arr[129] =
      (arrold[121] + 7 * arrold[123] - 77 * arrold[125] + 245 * arrold[127] - 358 * arrold[129] +
       182 * arrold[131]) * coeff;
      arr[130] =
      (4 * arrold[122] - 32 * arrold[124] + 108 * arrold[126] - 192 * arrold[128] +
       168 * arrold[130]) * coeff;
      arr[131] =
      (arrold[121] - 9 * arrold[123] + 35 * arrold[125] - 75 * arrold[127] + 90 * arrold[129] -
       42 * arrold[131]) * coeff;
      arr[132] =
      (-462 * arrold[133] + 330 * arrold[135] - 165 * arrold[137] + 55 * arrold[139] - 11 * arrold[141] +
       1 * arrold[143]) * coeff;
      arr[133] =
      (-252 * arrold[132] + 420 * arrold[134] - 240 * arrold[136] + 90 * arrold[138] - 20 * arrold[140] +
       2 * arrold[142]) * coeff;
      arr[134] =
      (546 * arrold[133] - 342 * arrold[135] + 123 * arrold[137] - 17 * arrold[139] - 3 * arrold[141] +
       1 * arrold[143]) * coeff;
      arr[135] =
      (364 * arrold[132] - 532 * arrold[134] + 176 * arrold[136] + 14 * arrold[138] - 28 * arrold[140] +
       6 * arrold[142]) * coeff;
      arr[136] =
      (-910 * arrold[133] + 330 * arrold[135] + 91 * arrold[137] - 105 * arrold[139] + 21 * arrold[141] +
       arrold[143]) * coeff;
      arr[137] =
      (-780 * arrold[132] + 820 * arrold[134] + 208 * arrold[136] - 318 * arrold[138] + 60 * arrold[140] +
       10 * arrold[142]) * coeff;
      arr[138] =
      (2210 * arrold[133] + 170 * arrold[135] - 901 * arrold[137] + 175 * arrold[139] + 61 * arrold[141] +
       arrold[143]) * coeff;
      arr[139] =
      (2652 * arrold[132] - 1156 * arrold[134] - 2448 * arrold[136] + 630 * arrold[138] + 308 * arrold[140] +
       14 * arrold[142]) * coeff;
      arr[140] =
      (-8398 * arrold[133] - 5814 * arrold[135] + 2907 * arrold[137] + 1463 * arrold[139] +
       117 * arrold[141] + arrold[143]) * coeff;
      arr[141] =
      (-16796 * arrold[132] - 6460 * arrold[134] + 15504 * arrold[136] + 6954 * arrold[138] +
       780 * arrold[140] + 18 * arrold[142]) * coeff;
      arr[142] =
      (58786 * arrold[133] + 87210 * arrold[135] + 33915 * arrold[137] + 4655 * arrold[139] +
       189 * arrold[141] + arrold[143]) * coeff;
      arr[143] =
      (352716 * arrold[132] + 497420 * arrold[134] + 170544 * arrold[136] + 26334 * arrold[138] +
       1540 * arrold[140] + 22 * arrold[142]) * coeff;
#endif

#if EXPANSION > 11
      coeff = 1 / 2048.0;
      arr[144] =
      (24 * arrold[145] + 2024 * arrold[147] + 42504 * arrold[149] + 346104 * arrold[151] +
       1307504 * arrold[153] + 2496144 * arrold[155]) * coeff;
      arr[145] =
      (arrold[144] + 230 * arrold[146] + 7084 * arrold[148] + 67298 * arrold[150] + 245157 * arrold[152] +
       326876 * arrold[154]) * coeff;
      arr[146] =
      (20 * arrold[145] + 1100 * arrold[147] + 13244 * arrold[149] + 47652 * arrold[151] +
       28424 * arrold[153] - 90440 * arrold[155]) * coeff;
      arr[147] =
      (arrold[144] + 150 * arrold[146] + 2604 * arrold[148] + 9842 * arrold[150] - 2907 * arrold[152] -
       34884 * arrold[154]) * coeff;
      arr[148] =
      (16 * arrold[145] + 496 * arrold[147] + 2224 * arrold[149] - 2736 * arrold[151] - 10336 * arrold[153] +
       10336 * arrold[155]) * coeff;
      arr[149] =
      (arrold[144] + 86 * arrold[146] + 556 * arrold[148] - 1102 * arrold[150] - 2907 * arrold[152] +
       6460 * arrold[154]) * coeff;
      arr[150] =
      (12 * arrold[145] + 148 * arrold[147] - 348 * arrold[149] - 900 * arrold[151] + 3128 * arrold[153] -
       2040 * arrold[155]) * coeff;
      arr[151] =
      (arrold[144] + 38 * arrold[146] - 84 * arrold[148] - 350 * arrold[150] + 1445 * arrold[152] -
       1700 * arrold[154]) * coeff;
      arr[152] =
      (8 * arrold[145] - 8 * arrold[147] - 168 * arrold[149] + 680 * arrold[151] - 1072 * arrold[153] +
       560 * arrold[155]) * coeff;
      arr[153] =
      (arrold[144] + 6 * arrold[146] - 84 * arrold[148] + 322 * arrold[150] - 603 * arrold[152] +
       540 * arrold[154]) * coeff;
      arr[154] =
      (4 * arrold[145] - 36 * arrold[147] + 140 * arrold[149] - 300 * arrold[151] + 360 * arrold[153] -
       168 * arrold[155]) * coeff;
      arr[155] =
      (arrold[144] - 10 * arrold[146] + 44 * arrold[148] - 110 * arrold[150] + 165 * arrold[152] -
       132 * arrold[154]) * coeff;
      arr[156] =
      (462 * arrold[156] - 792 * arrold[158] + 495 * arrold[160] - 220 * arrold[162] + 66 * arrold[164] -
       12 * arrold[166] + 1 * arrold[168]) * coeff;
      arr[157] =
      (-924 * arrold[157] + 660 * arrold[159] - 330 * arrold[161] + 110 * arrold[163] - 22 * arrold[165] +
       2 * arrold[167]) * coeff;
      arr[158] =
      (-546 * arrold[156] + 888 * arrold[158] - 465 * arrold[160] + 140 * arrold[162] - 14 * arrold[164] -
       4 * arrold[166] + 1 * arrold[168]) * coeff;
      arr[159] =
      (1260 * arrold[157] - 708 * arrold[159] + 162 * arrold[161] + 42 * arrold[163] - 34 * arrold[165] +
       6 * arrold[167]) * coeff;
      arr[160] =
      (910 * arrold[156] - 1240 * arrold[158] + 239 * arrold[160] + 196 * arrold[162] - 126 * arrold[164] +
       20 * arrold[166] + 1 * arrold[168]) * coeff;
      arr[161] =
      (-2380 * arrold[157] + 612 * arrold[159] + 526 * arrold[161] - 378 * arrold[163] + 50 * arrold[165] +
       10 * arrold[167]) * coeff;
      arr[162] =
      (-2210 * arrold[156] + 2040 * arrold[158] + 1071 * arrold[160] - 1076 * arrold[162] +
       114 * arrold[164] + 60 * arrold[166] + 1 * arrold[168]) * coeff;
      arr[163] =
      (6460 * arrold[157] + 1292 * arrold[159] - 3078 * arrold[161] + 322 * arrold[163] + 294 * arrold[165] +
       14 * arrold[167]) * coeff;
      arr[164] =
      (8398 * arrold[156] - 2584 * arrold[158] - 8721 * arrold[160] + 1444 * arrold[162] +
       1346 * arrold[164] + 116 * arrold[166] + 1 * arrold[168]) * coeff;
      arr[165] =
      (-27132 * arrold[157] - 21964 * arrold[159] + 8550 * arrold[161] + 6174 * arrold[163] +
       762 * arrold[165] + 18 * arrold[167]) * coeff;
      arr[166] =
      (-58786 * arrold[156] - 28424 * arrold[158] + 53295 * arrold[160] + 29260 * arrold[162] +
       4466 * arrold[164] + 188 * arrold[166] + 1 * arrold[168]) * coeff;
      arr[167] =
      (208012 * arrold[157] + 326876 * arrold[159] + 144210 * arrold[161] + 24794 * arrold[163] +
       1518 * arrold[165] + 22 * arrold[167]) * coeff;
      arr[168] =
      (1352078 * arrold[156] + 1961256 * arrold[158] + 735471 * arrold[160] + 134596 * arrold[162] +
       10626 * arrold[164] + 276 * arrold[166] + arrold[168]) * coeff;
#endif

#if EXPANSION > 12
      coeff = 1 / 4096.0;
      arr[169] =
      (arrold[169] + 325 * arrold[171] + 14950 * arrold[173] + 230230 * arrold[175] + 1562275 * arrold[177] +
       5311735 * arrold[179] + 9657700 * arrold[181]) * coeff;
      arr[170] =
      (24 * arrold[170] + 2000 * arrold[172] + 40480 * arrold[174] + 303600 * arrold[176] +
       961400 * arrold[178] + 1188640 * arrold[180]) * coeff;
      arr[171] =
      (arrold[169] + 229 * arrold[171] + 6854 * arrold[173] + 60214 * arrold[175] + 177859 * arrold[177] +
       81719 * arrold[179] - 326876 * arrold[181]) * coeff;
      arr[172] =
      (20 * arrold[170] + 1080 * arrold[172] + 12144 * arrold[174] + 34408 * arrold[176] -
       19228 * arrold[178] - 118864 * arrold[180]) * coeff;
      arr[173] =
      (arrold[169] + 149 * arrold[171] + 2454 * arrold[173] + 7238 * arrold[175] - 12749 * arrold[177] -
       31977 * arrold[179] + 34884 * arrold[181]) * coeff;
      arr[174] =
      (16 * arrold[170] + 480 * arrold[172] + 1728 * arrold[174] - 4960 * arrold[176] - 7600 * arrold[178] +
       20672 * arrold[180]) * coeff;
      arr[175] =
      (arrold[169] + 85 * arrold[171] + 470 * arrold[173] - 1658 * arrold[175] - 1805 * arrold[177] +
       9367 * arrold[179] - 6460 * arrold[181]) * coeff;
      arr[176] =
      (12 * arrold[170] + 136 * arrold[172] - 496 * arrold[174] - 552 * arrold[176] + 4028 * arrold[178] -
       5168 * arrold[180]) * coeff;
      arr[177] =
      (arrold[169] + 37 * arrold[171] - 122 * arrold[173] - 266 * arrold[175] + 1795 * arrold[177] -
       3145 * arrold[179] + 1700 * arrold[181]) * coeff;
      arr[178] =
      (8 * arrold[170] - 16 * arrold[172] - 160 * arrold[174] + 848 * arrold[176] - 1752 * arrold[178] +
       1632 * arrold[180]) * coeff;
      arr[179] =
      (arrold[169] + 5 * arrold[171] - 90 * arrold[173] + 406 * arrold[175] - 925 * arrold[177] +
       1143 * arrold[179] - 540 * arrold[181]) * coeff;
      arr[180] =
      (4 * arrold[170] - 40 * arrold[172] + 176 * arrold[174] - 440 * arrold[176] + 660 * arrold[178] -
       528 * arrold[180]) * coeff;
      arr[181] =
      (arrold[169] - 11 * arrold[171] + 54 * arrold[173] - 154 * arrold[175] + 275 * arrold[177] -
       297 * arrold[179] + 132 * arrold[181]) * coeff;
      arr[182] =
      (1716 * arrold[183] - 1287 * arrold[185] + 715 * arrold[187] - 286 * arrold[189] + 78 * arrold[191] -
       13 * arrold[193] + arrold[195]) * coeff;
      arr[183] =
      (924 * arrold[182] - 1584 * arrold[184] + 990 * arrold[186] - 440 * arrold[188] + 132 * arrold[190] -
       24 * arrold[192] + 2 * arrold[194]) * coeff;
      arr[184] =
      (-1980 * arrold[183] + 1353 * arrold[185] - 605 * arrold[187] + 154 * arrold[189] - 10 * arrold[191] -
       5 * arrold[193] + arrold[195]) * coeff;
      arr[185] =
      (-1260 * arrold[182] + 1968 * arrold[184] - 870 * arrold[186] + 120 * arrold[188] + 76 * arrold[190] -
       40 * arrold[192] + 6 * arrold[194]) * coeff;
      arr[186] =
      (3060 * arrold[183] - 1479 * arrold[185] + 43 * arrold[187] + 322 * arrold[189] - 146 * arrold[191] +
       19 * arrold[193] + arrold[195]) * coeff;
      arr[187] =
      (2380 * arrold[182] - 2992 * arrold[184] + 86 * arrold[186] + 904 * arrold[188] - 428 * arrold[190] +
       40 * arrold[192] + 10 * arrold[194]) * coeff;
      arr[188] =
      (-6460 * arrold[183] + 969 * arrold[185] + 2147 * arrold[187] - 1190 * arrold[189] + 54 * arrold[191] +
       59 * arrold[193] + arrold[195]) * coeff;
      arr[189] =
      (-6460 * arrold[182] + 5168 * arrold[184] + 4370 * arrold[186] - 3400 * arrold[188] + 28 * arrold[190] +
       280 * arrold[192] + 14 * arrold[194]) * coeff;
      arr[190] =
      (19380 * arrold[183] + 6137 * arrold[185] - 10165 * arrold[187] + 98 * arrold[189] +
       1230 * arrold[191] + 115 * arrold[193] + arrold[195]) * coeff;
      arr[191] =
      (27132 * arrold[182] - 5168 * arrold[184] - 30514 * arrold[186] + 2376 * arrold[188] +
       5412 * arrold[190] + 744 * arrold[192] + 18 * arrold[194]) * coeff;
      arr[192] =
      (-89148 * arrold[183] - 81719 * arrold[185] + 24035 * arrold[187] + 24794 * arrold[189] +
       4278 * arrold[191] + 187 * arrold[193] + arrold[195]) * coeff;
      arr[193] =
      (-208012 * arrold[182] - 118864 * arrold[184] + 182666 * arrold[186] + 119416 * arrold[188] +
       23276 * arrold[190] + 1496 * arrold[192] + 22 * arrold[194]) * coeff;
      arr[194] =
      (742900 * arrold[183] + 1225785 * arrold[185] + 600875 * arrold[187] + 123970 * arrold[189] +
       10350 * arrold[191] + 275 * arrold[193] + arrold[195]) * coeff;
      arr[195] =
      (5200300 * arrold[182] + 7726160 * arrold[184] + 3124550 * arrold[186] + 657800 * arrold[188] +
       65780 * arrold[190] + 2600 * arrold[192] + 26 * arrold[194]) * coeff;
#endif

#if EXPANSION > 13
      coeff = 1.0 / 8192.0;
      arr[196] =
      (28 * arrold[197] + 3276 * arrold[199] + 98280 * arrold[201] + 1184040 * arrold[203] +
       6906900 * arrold[205] + 21474180 * arrold[207] + 37442160 * arrold[209]) * coeff;
      arr[197] =
      (arrold[196] + 324 * arrold[198] + 14625 * arrold[200] + 215280 * arrold[202] + 1332045 * arrold[204] +
       3749460 * arrold[206] + 4345965 * arrold[208]) * coeff;
      arr[198] =
      (24 * arrold[197] + 1976 * arrold[199] + 38480 * arrold[201] + 263120 * arrold[203] +
       657800 * arrold[205] + 227240 * arrold[207] - 1188640 * arrold[209]) * coeff;
      arr[199] =
      (arrold[196] + 228 * arrold[198] + 6625 * arrold[200] + 53360 * arrold[202] + 117645 * arrold[204] -
       96140 * arrold[206] - 408595 * arrold[208]) * coeff;
      arr[200] =
      (20 * arrold[197] + 1060 * arrold[199] + 11064 * arrold[201] + 22264 * arrold[203] -
       53636 * arrold[205] - 99636 * arrold[207] + 118864 * arrold[209]) * coeff;
      arr[201] =
      (arrold[196] + 148 * arrold[198] + 2305 * arrold[200] + 4784 * arrold[202] - 19987 * arrold[204] -
       19228 * arrold[206] + 66861 * arrold[208]) * coeff;
      arr[202] =
      (16 * arrold[197] + 464 * arrold[199] + 1248 * arrold[201] - 6688 * arrold[203] - 2640 * arrold[205] +
       28272 * arrold[207] - 20672 * arrold[209]) * coeff;
      arr[203] =
      (arrold[196] + 84 * arrold[198] + 385 * arrold[200] - 2128 * arrold[202] - 147 * arrold[204] +
       11172 * arrold[206] - 15827 * arrold[208]) * coeff;
      arr[204] =
      (12 * arrold[197] + 124 * arrold[199] - 632 * arrold[201] - 56 * arrold[203] + 4580 * arrold[205] -
       9196 * arrold[207] + 5168 * arrold[209]) * coeff;
      arr[205] =
      (arrold[196] + 36 * arrold[198] - 159 * arrold[200] - 144 * arrold[202] + 2061 * arrold[204] -
       4940 * arrold[206] + 4845 * arrold[208]) * coeff;
      arr[206] =
      (8 * arrold[197] - 24 * arrold[199] - 144 * arrold[201] + 1008 * arrold[203] - 2600 * arrold[205] +
       3384 * arrold[207] - 1632 * arrold[209]) * coeff;
      arr[207] =
      (arrold[196] + 4 * arrold[198] - 95 * arrold[200] + 496 * arrold[202] - 1331 * arrold[204] +
       2068 * arrold[206] - 1683 * arrold[208]) * coeff;
      arr[208] =
      (4 * arrold[197] - 44 * arrold[199] + 216 * arrold[201] - 616 * arrold[203] + 1100 * arrold[205] -
       1188 * arrold[207] + 528 * arrold[209]) * coeff;
      arr[209] =
      (arrold[196] - 12 * arrold[198] + 65 * arrold[200] - 208 * arrold[202] + 429 * arrold[204] -
       572 * arrold[206] + 429 * arrold[208]) * coeff;
      arr[210] =
      (-1716 * arrold[210] + 3003 * arrold[212] - 2002 * arrold[214] + 1001 * arrold[216] -
       364 * arrold[218] + 91 * arrold[220] - 14 * arrold[222] + arrold[224]) * coeff;
      arr[211] =
      (3432 * arrold[211] - 2574 * arrold[213] + 1430 * arrold[215] - 572 * arrold[217] + 156 * arrold[219] -
       26 * arrold[221] + 2 * arrold[223]) * coeff;
      arr[212] =
      (1980 * arrold[210] - 3333 * arrold[212] + 1958 * arrold[214] - 759 * arrold[216] + 164 * arrold[218] -
       5 * arrold[220] - 6 * arrold[222] + arrold[224]) * coeff;
      arr[213] =
      (-4488 * arrold[211] + 2838 * arrold[213] - 990 * arrold[215] + 44 * arrold[217] + 116 * arrold[219] -
       46 * arrold[221] + 6 * arrold[223]) * coeff;
      arr[214] =
      (-3060 * arrold[210] + 4539 * arrold[212] - 1522 * arrold[214] - 279 * arrold[216] + 468 * arrold[218] -
       165 * arrold[220] + 18 * arrold[222] + arrold[224]) * coeff;
      arr[215] =
      (7752 * arrold[211] - 3078 * arrold[213] - 818 * arrold[215] + 1332 * arrold[217] - 468 * arrold[219] +
       30 * arrold[221] + 10 * arrold[223]) * coeff;
      arr[216] =
      (6460 * arrold[210] - 7429 * arrold[212] - 1178 * arrold[214] + 3337 * arrold[216] -
       1244 * arrold[218] - 5 * arrold[220] + 58 * arrold[222] + arrold[224]) * coeff;
      arr[217] =
      (-18088 * arrold[211] + 798 * arrold[213] + 7770 * arrold[215] - 3428 * arrold[217] -
       252 * arrold[219] + 266 * arrold[221] + 14 * arrold[223]) * coeff;
      arr[218] =
      (-19380 * arrold[210] + 13243 * arrold[212] + 16302 * arrold[214] - 10263 * arrold[216] -
       1132 * arrold[218] + 1115 * arrold[220] + 114 * arrold[222] + arrold[224]) * coeff;
      arr[219] =
      (59432 * arrold[211] + 25346 * arrold[213] - 32890 * arrold[215] - 3036 * arrold[217] +
       4668 * arrold[219] + 726 * arrold[221] + 18 * arrold[223]) * coeff;
      arr[220] =
      (89148 * arrold[210] - 7429 * arrold[212] - 105754 * arrold[214] - 759 * arrold[216] +
       20516 * arrold[218] + 4091 * arrold[220] + 186 * arrold[222] + arrold[224]) * coeff;
      arr[221] =
      (-297160 * arrold[211] - 301530 * arrold[213] + 63250 * arrold[215] + 96140 * arrold[217] +
       21780 * arrold[219] + 1474 * arrold[221] + 22 * arrold[223]) * coeff;
      arr[222] =
      (-742900 * arrold[210] - 482885 * arrold[212] + 624910 * arrold[214] + 476905 * arrold[216] +
       113620 * arrold[218] + 10075 * arrold[220] + 274 * arrold[222] + arrold[224]) * coeff;
      arr[223] =
      (2674440 * arrold[211] + 4601610 * arrold[213] + 2466750 * arrold[215] + 592020 * arrold[217] +
       63180 * arrold[219] + 2574 * arrold[221] + 26 * arrold[223]) * coeff;
      arr[224] =
      (20058300 * arrold[210] + 30421755 * arrold[212] + 13123110 * arrold[214] + 3108105 * arrold[216] +
       376740 * arrold[218] + 20475 * arrold[220] + 378 * arrold[222] + arrold[224]) * coeff;
#endif

#if EXPANSION > 14
      coeff = 1.0 / 16384.0;
      arr[225] =
      (arrold[225] + 435 * arrold[227] + 27405 * arrold[229] + 593775 * arrold[231] + 5852925 * arrold[233] +
       30045015 * arrold[235] + 86493225 * arrold[237] + 145422675 * arrold[239]) * coeff;
      arr[226] =
      (28 * arrold[226] + 3248 * arrold[228] + 95004 * arrold[230] + 1085760 * arrold[232] +
       5722860 * arrold[234] + 14567280 * arrold[236] + 15967980 * arrold[238]) * coeff;
      arr[227] =
      (arrold[225] + 323 * arrold[227] + 14301 * arrold[229] + 200655 * arrold[231] + 1116765 * arrold[233] +
       2417415 * arrold[235] + 596505 * arrold[237] - 4345965 * arrold[239]) * coeff;
      arr[228] =
      (24 * arrold[226] + 1952 * arrold[228] + 36504 * arrold[230] + 224640 * arrold[232] +
       394680 * arrold[234] - 430560 * arrold[236] - 1415880 * arrold[238]) * coeff;
      arr[229] =
      (arrold[225] + 227 * arrold[227] + 6397 * arrold[229] + 46735 * arrold[231] + 64285 * arrold[233] -
       213785 * arrold[235] - 312455 * arrold[237] + 408595 * arrold[239]) * coeff;
      arr[230] =
      (20 * arrold[226] + 1040 * arrold[228] + 10004 * arrold[230] + 11200 * arrold[232] -
       75900 * arrold[234] - 46000 * arrold[236] + 218500 * arrold[238]) * coeff;
      arr[231] =
      (arrold[225] + 147 * arrold[227] + 2157 * arrold[229] + 2479 * arrold[231] - 24771 * arrold[233] +
       759 * arrold[235] + 86089 * arrold[237] - 66861 * arrold[239]) * coeff;
      arr[232] =
      (16 * arrold[226] + 448 * arrold[228] + 784 * arrold[230] - 7936 * arrold[232] + 4048 * arrold[234] +
       30912 * arrold[236] - 48944 * arrold[238]) * coeff;
      arr[233] =
      (arrold[225] + 83 * arrold[227] + 301 * arrold[229] - 2513 * arrold[231] + 1981 * arrold[233] +
       11319 * arrold[235] - 26999 * arrold[237] + 15827 * arrold[239]) * coeff;
      arr[234] =
      (12 * arrold[226] + 112 * arrold[228] - 756 * arrold[230] + 576 * arrold[232] + 4636 * arrold[234] -
       13776 * arrold[236] + 14364 * arrold[238]) * coeff;
      arr[235] =
      (arrold[225] + 35 * arrold[227] - 195 * arrold[229] + 15 * arrold[231] + 2205 * arrold[233] -
       7001 * arrold[235] + 9785 * arrold[237] - 4845 * arrold[239]) * coeff;
      arr[236] =
      (8 * arrold[226] - 32 * arrold[228] - 120 * arrold[230] + 1152 * arrold[232] - 3608 * arrold[234] +
       5984 * arrold[236] - 5016 * arrold[238]) * coeff;
      arr[237] =
      (arrold[225] + 3 * arrold[227] - 99 * arrold[229] + 591 * arrold[231] - 1827 * arrold[233] +
       3399 * arrold[235] - 3751 * arrold[237] + 1683 * arrold[239]) * coeff;
      arr[238] =
      (4 * arrold[226] - 48 * arrold[228] + 260 * arrold[230] - 832 * arrold[232] + 1716 * arrold[234] -
       2288 * arrold[236] + 1716 * arrold[238]) * coeff;
      arr[239] =
      (arrold[225] - 13 * arrold[227] + 77 * arrold[229] - 273 * arrold[231] + 637 * arrold[233] -
       1001 * arrold[235] + 1001 * arrold[237] - 429 * arrold[239]) * coeff;
      arr[240] =
      (-6435 * arrold[241] + 5005 * arrold[243] - 3003 * arrold[245] + 1365 * arrold[247] -
       455 * arrold[249] + 105 * arrold[251] - 15 * arrold[253] + arrold[255]) * coeff;
      arr[241] =
      (-3432 * arrold[240] + 6006 * arrold[242] - 4004 * arrold[244] + 2002 * arrold[246] -
       728 * arrold[248] + 182 * arrold[250] - 28 * arrold[252] + 2 * arrold[254]) * coeff;
      arr[242] =
      (7293 * arrold[241] - 5291 * arrold[243] + 2717 * arrold[245] - 923 * arrold[247] + 169 * arrold[249] +
       arrold[251] - 7 * arrold[253] + arrold[255]) * coeff;
      arr[243] =
      (4488 * arrold[240] - 7326 * arrold[242] + 3828 * arrold[244] - 1034 * arrold[246] - 72 * arrold[248] +
       162 * arrold[250] - 52 * arrold[252] + 6 * arrold[254]) * coeff;
      arr[244] =
      (-10659 * arrold[241] + 6061 * arrold[243] - 1243 * arrold[245] - 747 * arrold[247] +
       633 * arrold[249] - 183 * arrold[251] + 17 * arrold[253] + arrold[255]) * coeff;
      arr[245] =
      (-7752 * arrold[240] + 10830 * arrold[242] - 2260 * arrold[244] - 2150 * arrold[246] +
       1800 * arrold[248] - 498 * arrold[250] + 20 * arrold[252] + 10 * arrold[254]) * coeff;
      arr[246] =
      (20349 * arrold[241] - 6251 * arrold[243] - 4515 * arrold[245] + 4581 * arrold[247] -
       1239 * arrold[249] - 63 * arrold[251] + 57 * arrold[253] + arrold[255]) * coeff;
      arr[247] =
      (18088 * arrold[240] - 18886 * arrold[242] - 6972 * arrold[244] + 11198 * arrold[246] -
       3176 * arrold[248] - 518 * arrold[250] + 252 * arrold[252] + 14 * arrold[254]) * coeff;
      arr[248] =
      (-52003 * arrold[241] - 3059 * arrold[243] + 26565 * arrold[245] - 9131 * arrold[247] -
       2247 * arrold[249] + 1001 * arrold[251] + 113 * arrold[253] + arrold[255]) * coeff;
      arr[249] =
      (-59432 * arrold[240] + 34086 * arrold[242] + 58236 * arrold[244] - 29854 * arrold[246] -
       7704 * arrold[248] + 3942 * arrold[250] + 708 * arrold[252] + 18 * arrold[254]) * coeff;
      arr[250] =
      (185725 * arrold[241] + 98325 * arrold[243] - 104995 * arrold[245] - 21275 * arrold[247] +
       16425 * arrold[249] + 3905 * arrold[251] + 185 * arrold[253] + arrold[255]) * coeff;
      arr[251] =
      (297160 * arrold[240] + 4370 * arrold[242] - 364780 * arrold[244] - 32890 * arrold[246] +
       74360 * arrold[248] + 20306 * arrold[250] + 1452 * arrold[252] + 22 * arrold[254]) * coeff;
      arr[252] =
      (-1002915 * arrold[241] - 1107795 * arrold[243] + 148005 * arrold[245] + 363285 * arrold[247] +
       103545 * arrold[249] + 9801 * arrold[251] + 273 * arrold[253] + arrold[255]) * coeff;
      arr[253] =
      (-2674440 * arrold[240] - 1927170 * arrold[242] + 2134860 * arrold[244] + 1874730 * arrold[246] +
       528840 * arrold[248] + 60606 * arrold[250] + 2548 * arrold[252] + 26 * arrold[254]) * coeff;
      arr[254] =
      (9694845 * arrold[241] + 17298645 * arrold[243] + 10015005 * arrold[245] + 2731365 * arrold[247] +
       356265 * arrold[249] + 20097 * arrold[251] + 377 * arrold[253] + arrold[255]) * coeff;
      arr[255] =
      (77558760 * arrold[240] + 119759850 * arrold[242] + 54627300 * arrold[244] + 14307150 * arrold[246] +
       2035800 * arrold[248] + 142506 * arrold[250] + 4060 * arrold[252] + 30 * arrold[254]) * coeff;
#endif

#if EXPANSION > 15
      coeff = 1.0 / 32768.0;
      arr[256] =
      (32 * arrold[257] + 4960 * arrold[259] + 201376 * arrold[261] + 3365856 * arrold[263] +
       28048800 * arrold[265] + 129024480 * arrold[267] + 347373600 * arrold[269] +
       565722720 * arrold[271]) * coeff;
      arr[257] =
      (arrold[256] + 434 * arrold[258] + 26970 * arrold[260] + 566370 * arrold[262] + 5259150 * arrold[264] +
       24192090 * arrold[266] + 56448210 * arrold[268] + 58929450 * arrold[270]) * coeff;
      arr[258] =
      (28 * arrold[257] + 3220 * arrold[259] + 91756 * arrold[261] + 990756 * arrold[263] +
       4637100 * arrold[265] + 8844420 * arrold[267] + 1400700 * arrold[269] -
       15967980 * arrold[271]) * coeff;
      arr[259] =
      (arrold[256] + 322 * arrold[258] + 13978 * arrold[260] + 186354 * arrold[262] + 916110 * arrold[264] +
       1300650 * arrold[266] - 1820910 * arrold[268] - 4942470 * arrold[270]) * coeff;
      arr[260] =
      (24 * arrold[257] + 1928 * arrold[259] + 34552 * arrold[261] + 188136 * arrold[263] +
       170040 * arrold[265] - 825240 * arrold[267] - 985320 * arrold[269] + 1415880 * arrold[271]) * coeff;
      arr[261] =
      (arrold[256] + 226 * arrold[258] + 6170 * arrold[260] + 40338 * arrold[262] + 17550 * arrold[264] -
       278070 * arrold[266] - 98670 * arrold[268] + 721050 * arrold[270]) * coeff;
      arr[262] =
      (20 * arrold[257] + 1020 * arrold[259] + 8964 * arrold[261] + 1196 * arrold[263] - 87100 * arrold[265] +
       29900 * arrold[267] + 264500 * arrold[269] - 218500 * arrold[271]) * coeff;
      arr[263] =
      (arrold[256] + 146 * arrold[258] + 2010 * arrold[260] + 322 * arrold[262] - 27250 * arrold[264] +
       25530 * arrold[266] + 85330 * arrold[268] - 152950 * arrold[270]) * coeff;
      arr[264] =
      (16 * arrold[257] + 432 * arrold[259] + 336 * arrold[261] - 8720 * arrold[263] + 11984 * arrold[265] +
       26864 * arrold[267] - 79856 * arrold[269] + 48944 * arrold[271]) * coeff;
      arr[265] =
      (arrold[256] + 82 * arrold[258] + 218 * arrold[260] - 2814 * arrold[262] + 4494 * arrold[264] +
       9338 * arrold[266] - 38318 * arrold[268] + 42826 * arrold[270]) * coeff;
      arr[266] =
      (12 * arrold[257] + 100 * arrold[259] - 868 * arrold[261] + 1332 * arrold[263] + 4060 * arrold[265] -
       18412 * arrold[267] + 28140 * arrold[269] - 14364 * arrold[271]) * coeff;
      arr[267] =
      (arrold[256] + 34 * arrold[258] - 230 * arrold[260] + 210 * arrold[262] + 2190 * arrold[264] -
       9206 * arrold[266] + 16786 * arrold[268] - 14630 * arrold[270]) * coeff;
      arr[268] =
      (8 * arrold[257] - 40 * arrold[259] - 88 * arrold[261] + 1272 * arrold[263] - 4760 * arrold[265] +
       9592 * arrold[267] - 11000 * arrold[269] + 5016 * arrold[271]) * coeff;
      arr[269] =
      (arrold[256] + 2 * arrold[258] - 102 * arrold[260] + 690 * arrold[262] - 2418 * arrold[264] +
       5226 * arrold[266] - 7150 * arrold[268] + 5434 * arrold[270]) * coeff;
      arr[270] =
      (4 * arrold[257] - 52 * arrold[259] + 308 * arrold[261] - 1092 * arrold[263] + 2548 * arrold[265] -
       4004 * arrold[267] + 4004 * arrold[269] - 1716 * arrold[271]) * coeff;
      arr[271] =
      (arrold[256] - 14 * arrold[258] + 90 * arrold[260] - 350 * arrold[262] + 910 * arrold[264] -
       1638 * arrold[266] + 2002 * arrold[268] - 1430 * arrold[270]) * coeff;
      arr[272] =
      (6435 * arrold[272] - 11440 * arrold[274] + 8008 * arrold[276] - 4368 * arrold[278] +
       1820 * arrold[280] - 560 * arrold[282] + 120 * arrold[284] - 16 * arrold[286] + arrold[288]) * coeff;
      arr[273] =
      (-12870 * arrold[273] + 10010 * arrold[275] - 6006 * arrold[277] + 2730 * arrold[279] -
       910 * arrold[281] + 210 * arrold[283] - 30 * arrold[285] + 2 * arrold[287]) * coeff;
      arr[274] =
      (-7293 * arrold[272] + 12584 * arrold[274] - 8008 * arrold[276] + 3640 * arrold[278] -
       1092 * arrold[280] + 168 * arrold[282] + 8 * arrold[284] - 8 * arrold[286] + arrold[288]) * coeff;
      arr[275] =
      (16302 * arrold[273] - 11154 * arrold[275] + 4862 * arrold[277] - 962 * arrold[279] -
       234 * arrold[281] + 214 * arrold[283] - 58 * arrold[285] + 6 * arrold[287]) * coeff;
      arr[276] =
      (10659 * arrold[272] - 16720 * arrold[274] + 7304 * arrold[276] - 496 * arrold[278] -
       1380 * arrold[280] + 816 * arrold[282] - 200 * arrold[284] + 16 * arrold[286] + arrold[288]) * coeff;
      arr[277] =
      (-26334 * arrold[273] + 13090 * arrold[275] - 110 * arrold[277] - 3950 * arrold[279] +
       2298 * arrold[281] - 518 * arrold[283] + 10 * arrold[285] + 10 * arrold[287]) * coeff;
      arr[278] =
      (-20349 * arrold[272] + 26600 * arrold[274] - 1736 * arrold[276] - 9096 * arrold[278] +
       5820 * arrold[280] - 1176 * arrold[282] - 120 * arrold[284] + 56 * arrold[286] + arrold[288]) * coeff;
      arr[279] =
      (55062 * arrold[273] - 11914 * arrold[275] - 18170 * arrold[277] + 14374 * arrold[279] -
       2658 * arrold[281] - 770 * arrold[283] + 238 * arrold[285] + 14 * arrold[287]) * coeff;
      arr[280] =
      (52003 * arrold[272] - 48944 * arrold[274] - 29624 * arrold[276] + 35696 * arrold[278] -
       6884 * arrold[280] - 3248 * arrold[282] + 888 * arrold[284] + 112 * arrold[286] + arrold[288]) * coeff;
      arr[281] =
      (-152950 * arrold[273] - 24150 * arrold[275] + 88090 * arrold[277] - 22150 * arrold[279] -
       11646 * arrold[281] + 3234 * arrold[283] + 690 * arrold[285] + 18 * arrold[287]) * coeff;
      arr[282] =
      (-185725 * arrold[272] + 87400 * arrold[274] + 203320 * arrold[276] - 83720 * arrold[278] -
       37700 * arrold[280] + 12520 * arrold[282] + 3720 * arrold[284] + 184 * arrold[286] +
       arrold[288]) * coeff;
      arr[283] =
      (589950 * arrold[273] + 369150 * arrold[275] - 331890 * arrold[277] - 107250 * arrold[279] +
       54054 * arrold[281] + 18854 * arrold[283] + 1430 * arrold[285] + 22 * arrold[287]) * coeff;
      arr[284] =
      (1002915 * arrold[272] + 104880 * arrold[274] - 1255800 * arrold[276] - 215280 * arrold[278] +
       259740 * arrold[280] + 93744 * arrold[282] + 9528 * arrold[284] + 272 * arrold[286] +
       arrold[288]) * coeff;
      arr[285] =
      (-3421710 * arrold[273] - 4062030 * arrold[275] + 260130 * arrold[277] + 1345890 * arrold[279] +
       468234 * arrold[281] + 58058 * arrold[283] + 2522 * arrold[285] + 26 * arrold[287]) * coeff;
      arr[286] =
      (-9694845 * arrold[272] - 7603800 * arrold[274] + 7283640 * arrold[276] + 7283640 * arrold[278] +
       2375100 * arrold[280] + 336168 * arrold[282] + 19720 * arrold[284] + 376 * arrold[286] +
       arrold[288]) * coeff;
      arr[287] =
      (35357670 * arrold[273] + 65132550 * arrold[275] + 40320150 * arrold[277] + 12271350 * arrold[279] +
       1893294 * arrold[281] + 138446 * arrold[283] + 4030 * arrold[285] + 30 * arrold[287]) * coeff;
      arr[288] =
      (300540195 * arrold[272] + 471435600 * arrold[274] + 225792840 * arrold[276] + 64512240 * arrold[278] +
       10518300 * arrold[280] + 906192 * arrold[282] + 35960 * arrold[284] + 496 * arrold[286] +
       arrold[288]) * coeff;
#endif


#if EXPANSION > 16
      coeff = 1.0 / 65536.0;
      arr[289] =
      (arrold[289] + 561 * arrold[291] + 46376 * arrold[293] + 1344904 * arrold[295] +
       18156204 * arrold[297] + 131128140 * arrold[299] + 548354040 * arrold[301] + 1391975640 * arrold[303] +
       2203961430 * arrold[305]) * coeff;
      arr[290] =
      (32 * arrold[290] + 4928 * arrold[292] + 196416 * arrold[294] + 3164480 * arrold[296] +
       24682944 * arrold[298] + 100975680 * arrold[300] + 218349120 * arrold[302] +
       218349120 * arrold[304]) * coeff;
      arr[291] =
      (arrold[289] + 433 * arrold[291] + 26536 * arrold[293] + 539400 * arrold[295] + 4692780 * arrold[297] +
       18932940 * arrold[299] + 32256120 * arrold[301] + 2481240 * arrold[303] -
       58929450 * arrold[305]) * coeff;
      arr[292] =
      (28 * arrold[290] + 3192 * arrold[292] + 88536 * arrold[294] + 899000 * arrold[296] +
       3646344 * arrold[298] + 4207320 * arrold[300] - 7443720 * arrold[302] -
       17368680 * arrold[304]) * coeff;
      arr[293] =
      (arrold[289] + 321 * arrold[291] + 13656 * arrold[293] + 172376 * arrold[295] + 729756 * arrold[297] +
       384540 * arrold[299] - 3121560 * arrold[301] - 3121560 * arrold[303] + 4942470 * arrold[305]) * coeff;
      arr[294] =
      (24 * arrold[290] + 1904 * arrold[292] + 32624 * arrold[294] + 153584 * arrold[296] -
       18096 * arrold[298] - 995280 * arrold[300] - 160080 * arrold[302] + 2401200 * arrold[304]) * coeff;
      arr[295] =
      (arrold[289] + 225 * arrold[291] + 5944 * arrold[293] + 34168 * arrold[295] - 22788 * arrold[297] -
       295620 * arrold[299] + 179400 * arrold[301] + 819720 * arrold[303] - 721050 * arrold[305]) * coeff;
      arr[296] =
      (20 * arrold[290] + 1000 * arrold[292] + 7944 * arrold[294] - 7768 * arrold[296] - 88296 * arrold[298] +
       117000 * arrold[300] + 234600 * arrold[302] - 483000 * arrold[304]) * coeff;
      arr[297] =
      (arrold[289] + 145 * arrold[291] + 1864 * arrold[293] - 1688 * arrold[295] - 27572 * arrold[297] +
       52780 * arrold[299] + 59800 * arrold[301] - 238280 * arrold[303] + 152950 * arrold[305]) * coeff;
      arr[298] =
      (16 * arrold[290] + 416 * arrold[292] - 96 * arrold[294] - 9056 * arrold[296] + 20704 * arrold[298] +
       14880 * arrold[300] - 106720 * arrold[302] + 128800 * arrold[304]) * coeff;
      arr[299] =
      (arrold[289] + 81 * arrold[291] + 136 * arrold[293] - 3032 * arrold[295] + 7308 * arrold[297] +
       4844 * arrold[299] - 47656 * arrold[301] + 81144 * arrold[303] - 42826 * arrold[305]) * coeff;
      arr[300] =
      (12 * arrold[290] + 88 * arrold[292] - 968 * arrold[294] + 2200 * arrold[296] + 2728 * arrold[298] -
       22472 * arrold[300] + 46552 * arrold[302] - 42504 * arrold[304]) * coeff;
      arr[301] =
      (arrold[289] + 33 * arrold[291] - 264 * arrold[293] + 440 * arrold[295] + 1980 * arrold[297] -
       11396 * arrold[299] + 25992 * arrold[301] - 31416 * arrold[303] + 14630 * arrold[305]) * coeff;
      arr[302] =
      (8 * arrold[290] - 48 * arrold[292] - 48 * arrold[294] + 1360 * arrold[296] - 6032 * arrold[298] +
       14352 * arrold[300] - 20592 * arrold[302] + 16016 * arrold[304]) * coeff;
      arr[303] =
      (arrold[289] + arrold[291] - 104 * arrold[293] + 792 * arrold[295] - 3108 * arrold[297] +
       7644 * arrold[299] - 12376 * arrold[301] + 12584 * arrold[303] - 5434 * arrold[305]) * coeff;
      arr[304] =
      (4 * arrold[290] - 56 * arrold[292] + 360 * arrold[294] - 1400 * arrold[296] + 3640 * arrold[298] -
       6552 * arrold[300] + 8008 * arrold[302] - 5720 * arrold[304]) * coeff;
      arr[305] =
      (arrold[289] - 15 * arrold[291] + 104 * arrold[293] - 440 * arrold[295] + 1260 * arrold[297] -
       2548 * arrold[299] + 3640 * arrold[301] - 3432 * arrold[303] + 1430 * arrold[305]) * coeff;
      arr[306] =
      (24310 * arrold[307] - 19448 * arrold[309] + 12376 * arrold[311] - 6188 * arrold[313] +
       2380 * arrold[315] - 680 * arrold[317] + 136 * arrold[319] - 17 * arrold[321] + arrold[323]) * coeff;
      arr[307] =
      (12870 * arrold[306] - 22880 * arrold[308] + 16016 * arrold[310] - 8736 * arrold[312] +
       3640 * arrold[314] - 1120 * arrold[316] + 240 * arrold[318] - 32 * arrold[320] +
       2 * arrold[322]) * coeff;
      arr[308] =
      (-27170 * arrold[307] + 20592 * arrold[309] - 11648 * arrold[311] + 4732 * arrold[313] -
       1260 * arrold[315] + 160 * arrold[317] + 16 * arrold[319] - 9 * arrold[321] + arrold[323]) * coeff;
      arr[309] =
      (-16302 * arrold[306] + 27456 * arrold[308] - 16016 * arrold[310] + 5824 * arrold[312] -
       728 * arrold[314] - 448 * arrold[316] + 272 * arrold[318] - 64 * arrold[320] +
       6 * arrold[322]) * coeff;
      arr[310] =
      (38038 * arrold[307] - 24024 * arrold[309] + 7800 * arrold[311] + 884 * arrold[313] -
       2196 * arrold[315] + 1016 * arrold[317] - 216 * arrold[319] + 15 * arrold[321] + arrold[323]) * coeff;
      arr[311] =
      (26334 * arrold[306] - 39424 * arrold[308] + 13200 * arrold[310] + 3840 * arrold[312] -
       6248 * arrold[314] + 2816 * arrold[316] - 528 * arrold[318] + 10 * arrold[322]) * coeff;
      arr[312] =
      (-67298 * arrold[307] + 28336 * arrold[309] + 7360 * arrold[311] - 14916 * arrold[313] +
       6996 * arrold[315] - 1056 * arrold[317] - 176 * arrold[319] + 55 * arrold[321] + arrold[323]) * coeff;
      arr[313] =
      (-55062 * arrold[306] + 66976 * arrold[308] + 6256 * arrold[310] - 32544 * arrold[312] +
       17032 * arrold[314] - 1888 * arrold[316] - 1008 * arrold[318] + 224 * arrold[320] +
       14 * arrold[322]) * coeff;
      arr[314] =
      (152950 * arrold[307] - 19320 * arrold[309] - 65320 * arrold[311] + 42580 * arrold[313] -
       3636 * arrold[315] - 4136 * arrold[317] + 776 * arrold[319] + 111 * arrold[321] + arrold[323]) * coeff;
      arr[315] =
      (152950 * arrold[306] - 128800 * arrold[308] - 112240 * arrold[310] + 110240 * arrold[312] -
       10504 * arrold[314] - 14880 * arrold[316] + 2544 * arrold[318] + 672 * arrold[320] +
       18 * arrold[322]) * coeff;
      arr[316] =
      (-458850 * arrold[307] - 115920 * arrold[309] + 287040 * arrold[311] - 46020 * arrold[313] -
       50220 * arrold[315] + 8800 * arrold[317] + 3536 * arrold[319] + 183 * arrold[321] +
       arrold[323]) * coeff;
      arr[317] =
      (-589950 * arrold[306] + 220800 * arrold[308] + 701040 * arrold[310] - 224640 * arrold[312] -
       161304 * arrold[314] + 35200 * arrold[316] + 17424 * arrold[318] + 1408 * arrold[320] +
       22 * arrold[322]) * coeff;
      arr[318] =
      (1900950 * arrold[307] + 1360680 * arrold[309] - 1040520 * arrold[311] - 475020 * arrold[313] +
       165996 * arrold[315] + 84216 * arrold[317] + 9256 * arrold[319] + 271 * arrold[321] +
       arrold[323]) * coeff;
      arr[319] =
      (3421710 * arrold[306] + 640320 * arrold[308] - 4322160 * arrold[310] - 1085760 * arrold[312] +
       877656 * arrold[314] + 410176 * arrold[316] + 55536 * arrold[318] + 2496 * arrold[320] +
       26 * arrold[322]) * coeff;
      arr[320] =
      (-11785890 * arrold[307] - 14887440 * arrold[309] + 4908540 * arrold[313] + 2038932 * arrold[315] +
       316448 * arrold[317] + 19344 * arrold[319] + 375 * arrold[321] + arrold[323]) * coeff;
      arr[321] =
      (-35357670 * arrold[306] - 29774880 * arrold[308] + 24812400 * arrold[310] + 28048800 * arrold[312] +
       10378056 * arrold[314] + 1754848 * arrold[316] + 134416 * arrold[318] + 4000 * arrold[320] +
       30 * arrold[322]) * coeff;
      arr[322] =
      (129644790 * arrold[307] + 245642760 * arrold[309] + 161280600 * arrold[311] + 53993940 * arrold[313] +
       9612108 * arrold[315] + 870232 * arrold[317] + 35464 * arrold[319] + 495 * arrold[321] +
       arrold[323]) * coeff;
      arr[323] =
      (1166803110 * arrold[306] + 1855967520 * arrold[308] + 927983760 * arrold[310] +
       286097760 * arrold[312] + 52451256 * arrold[314] + 5379616 * arrold[316] + 278256 * arrold[318] +
       5984 * arrold[320] + 34 * arrold[322]) * coeff;
#endif

#if EXPANSION > 17
      coeff = 1.0 / 131072.0;
      arr[324] =
      (36 * arrold[325] + 7140 * arrold[327] + 376992 * arrold[329] + 8347680 * arrold[331] +
       94143280 * arrold[333] + 600805296 * arrold[335] + 2310789600 * arrold[337] +
       5567902560 * arrold[339] + 8597496600 * arrold[341]) * coeff;
      arr[325] =
      (arrold[324] + 560 * arrold[326] + 45815 * arrold[328] + 1298528 * arrold[330] +
       16811300 * arrold[332] + 112971936 * arrold[334] + 417225900 * arrold[336] + 843621600 * arrold[338] +
       811985790 * arrold[340]) * coeff;
      arr[326] =
      (32 * arrold[325] + 4896 * arrold[327] + 191488 * arrold[329] + 2968064 * arrold[331] +
       21518464 * arrold[333] + 76292736 * arrold[335] + 117373440 * arrold[337] -
       218349120 * arrold[341]) * coeff;
      arr[327] =
      (arrold[324] + 432 * arrold[326] + 26103 * arrold[328] + 512864 * arrold[330] + 4153380 * arrold[332] +
       14240160 * arrold[334] + 13323180 * arrold[336] - 29774880 * arrold[338] -
       61410690 * arrold[340]) * coeff;
      arr[328] =
      (28 * arrold[325] + 3164 * arrold[327] + 85344 * arrold[329] + 810464 * arrold[331] +
       2747344 * arrold[333] + 560976 * arrold[335] - 11651040 * arrold[337] - 9924960 * arrold[339] +
       17368680 * arrold[341]) * coeff;
      arr[329] =
      (arrold[324] + 320 * arrold[326] + 13335 * arrold[328] + 158720 * arrold[330] + 557380 * arrold[332] -
       345216 * arrold[334] - 3506100 * arrold[336] + 8064030 * arrold[340]) * coeff;
      arr[330] =
      (24 * arrold[325] + 1880 * arrold[327] + 30720 * arrold[329] + 120960 * arrold[331] -
       171680 * arrold[333] - 977184 * arrold[335] + 835200 * arrold[337] + 2561280 * arrold[339] -
       2401200 * arrold[341]) * coeff;
      arr[331] =
      (arrold[324] + 224 * arrold[326] + 5719 * arrold[328] + 28224 * arrold[330] - 56956 * arrold[332] -
       272832 * arrold[334] + 475020 * arrold[336] + 640320 * arrold[338] - 1540770 * arrold[340]) * coeff;
      arr[332] =
      (20 * arrold[325] + 980 * arrold[327] + 6944 * arrold[329] - 15712 * arrold[331] - 80528 * arrold[333] +
       205296 * arrold[335] + 117600 * arrold[337] - 717600 * arrold[339] + 483000 * arrold[341]) * coeff;
      arr[333] =
      (arrold[324] + 144 * arrold[326] + 1719 * arrold[328] - 3552 * arrold[330] - 25884 * arrold[332] +
       80352 * arrold[334] + 7020 * arrold[336] - 298080 * arrold[338] + 391230 * arrold[340]) * coeff;
      arr[334] =
      (16 * arrold[325] + 400 * arrold[327] - 512 * arrold[329] - 8960 * arrold[331] + 29760 * arrold[333] -
       5824 * arrold[335] - 121600 * arrold[337] + 235520 * arrold[339] - 128800 * arrold[341]) * coeff;
      arr[335] =
      (arrold[324] + 80 * arrold[326] + 55 * arrold[328] - 3168 * arrold[330] + 10340 * arrold[332] -
       2464 * arrold[334] - 52500 * arrold[336] + 128800 * arrold[338] - 123970 * arrold[340]) * coeff;
      arr[336] =
      (12 * arrold[325] + 76 * arrold[327] - 1056 * arrold[329] + 3168 * arrold[331] + 528 * arrold[333] -
       25200 * arrold[335] + 69024 * arrold[337] - 89056 * arrold[339] + 42504 * arrold[341]) * coeff;
      arr[337] =
      (arrold[324] + 32 * arrold[326] - 297 * arrold[328] + 704 * arrold[330] + 1540 * arrold[332] -
       13376 * arrold[334] + 37388 * arrold[336] - 57408 * arrold[338] + 46046 * arrold[340]) * coeff;
      arr[338] =
      (8 * arrold[325] - 56 * arrold[327] + 1408 * arrold[331] - 7392 * arrold[333] + 20384 * arrold[335] -
       34944 * arrold[337] + 36608 * arrold[339] - 16016 * arrold[341]) * coeff;
      arr[339] =
      (arrold[324] - 105 * arrold[328] + 896 * arrold[330] - 3900 * arrold[332] + 10752 * arrold[334] -
       20020 * arrold[336] + 24960 * arrold[338] - 18018 * arrold[340]) * coeff;
      arr[340] =
      (4 * arrold[325] - 60 * arrold[327] + 416 * arrold[329] - 1760 * arrold[331] + 5040 * arrold[333] -
       10192 * arrold[335] + 14560 * arrold[337] - 13728 * arrold[339] + 5720 * arrold[341]) * coeff;
      arr[341] =
      (arrold[324] - 16 * arrold[326] + 119 * arrold[328] - 544 * arrold[330] + 1700 * arrold[332] -
       3808 * arrold[334] + 6188 * arrold[336] - 7072 * arrold[338] + 4862 * arrold[340]) * coeff;
      arr[342] =
      (-24310 * arrold[342] + 43758 * arrold[344] - 31824 * arrold[346] + 18564 * arrold[348] -
       8568 * arrold[350] + 3060 * arrold[352] - 816 * arrold[354] + 153 * arrold[356] - 18 * arrold[358] +
       arrold[360]) * coeff;
      arr[343] =
      (48620 * arrold[343] - 38896 * arrold[345] + 24752 * arrold[347] - 12376 * arrold[349] +
       4760 * arrold[351] - 1360 * arrold[353] + 272 * arrold[355] - 34 * arrold[357] +
       2 * arrold[359]) * coeff;
      arr[344] =
      (27170 * arrold[342] - 47762 * arrold[344] + 32240 * arrold[346] - 16380 * arrold[348] +
       5992 * arrold[350] - 1420 * arrold[352] + 144 * arrold[354] + 25 * arrold[356] - 10 * arrold[358] +
       arrold[360]) * coeff;
      arr[345] =
      (-60060 * arrold[343] + 43472 * arrold[345] - 21840 * arrold[347] + 6552 * arrold[349] -
       280 * arrold[351] - 720 * arrold[353] + 336 * arrold[355] - 70 * arrold[357] +
       6 * arrold[359]) * coeff;
      arr[346] =
      (-38038 * arrold[342] + 62062 * arrold[344] - 31824 * arrold[346] + 6916 * arrold[348] +
       3080 * arrold[350] - 3212 * arrold[352] + 1232 * arrold[354] - 231 * arrold[356] + 14 * arrold[358] +
       arrold[360]) * coeff;
      arr[347] =
      (92092 * arrold[343] - 52624 * arrold[345] + 9360 * arrold[347] + 10088 * arrold[349] -
       9064 * arrold[351] + 3344 * arrold[353] - 528 * arrold[355] - 10 * arrold[357] +
       10 * arrold[359]) * coeff;
      arr[348] =
      (67298 * arrold[342] - 95634 * arrold[344] + 20976 * arrold[346] + 22276 * arrold[348] -
       21912 * arrold[350] + 8052 * arrold[352] - 880 * arrold[354] - 231 * arrold[356] + 54 * arrold[358] +
       arrold[360]) * coeff;
      arr[349] =
      (-177100 * arrold[343] + 60720 * arrold[345] + 38800 * arrold[347] - 49576 * arrold[349] +
       18920 * arrold[351] - 880 * arrold[353] - 1232 * arrold[355] + 210 * arrold[357] +
       14 * arrold[359]) * coeff;
      arr[350] =
      (-152950 * arrold[342] + 172270 * arrold[344] + 46000 * arrold[346] - 107900 * arrold[348] +
       46216 * arrold[350] + 500 * arrold[352] - 4912 * arrold[354] + 665 * arrold[356] + 110 * arrold[358] +
       arrold[360]) * coeff;
      arr[351] =
      (434700 * arrold[343] - 16560 * arrold[345] - 222480 * arrold[347] + 120744 * arrold[349] +
       4376 * arrold[351] - 17424 * arrold[353] + 1872 * arrold[355] + 654 * arrold[357] +
       18 * arrold[359]) * coeff;
      arr[352] =
      (458850 * arrold[342] - 342930 * arrold[344] - 402960 * arrold[346] + 333060 * arrold[348] +
       4200 * arrold[350] - 59020 * arrold[352] + 5264 * arrold[354] + 3353 * arrold[356] +
       182 * arrold[358] + arrold[360]) * coeff;
      arr[353] =
      (-1400700 * arrold[343] - 480240 * arrold[345] + 925680 * arrold[347] - 63336 * arrold[349] -
       196504 * arrold[351] + 17776 * arrold[353] + 16016 * arrold[355] + 1386 * arrold[357] +
       22 * arrold[359]) * coeff;
      arr[354] =
      (-1900950 * arrold[342] + 540270 * arrold[344] + 2401200 * arrold[346] - 565500 * arrold[348] -
       641016 * arrold[350] + 81780 * arrold[352] + 74960 * arrold[354] + 8985 * arrold[356] +
       270 * arrold[358] + arrold[360]) * coeff;
      arr[355] =
      (6203100 * arrold[343] + 4962480 * arrold[345] - 3236400 * arrold[347] - 1963416 * arrold[349] +
       467480 * arrold[351] + 354640 * arrold[353] + 53040 * arrold[355] + 2470 * arrold[357] +
       26 * arrold[359]) * coeff;
      arr[356] =
      (11785890 * arrold[342] + 3101550 * arrold[344] - 14887440 * arrold[346] - 4908540 * arrold[348] +
       2869608 * arrold[350] + 1722484 * arrold[352] + 297104 * arrold[354] + 18969 * arrold[356] +
       374 * arrold[358] + arrold[360]) * coeff;
      arr[357] =
      (-40940460 * arrold[343] - 54587280 * arrold[345] - 3236400 * arrold[347] + 17670744 * arrold[349] +
       8623208 * arrold[351] + 1620432 * arrold[353] + 130416 * arrold[355] + 3970 * arrold[357] +
       30 * arrold[359]) * coeff;
      arr[358] =
      (-129644790 * arrold[342] - 115997970 * arrold[344] + 84362160 * arrold[346] + 107286660 * arrold[348] +
       44381832 * arrold[350] + 8741876 * arrold[352] + 834768 * arrold[354] + 34969 * arrold[356] +
       494 * arrold[358] + arrold[360]) * coeff;
      arr[359] =
      (477638700 * arrold[343] + 927983760 * arrold[345] + 641886000 * arrold[347] + 233646504 * arrold[349] +
       47071640 * arrold[351] + 5101360 * arrold[353] + 272272 * arrold[355] + 5950 * arrold[357] +
       34 * arrold[359]) * coeff;
      arr[360] =
      (4537567650 * arrold[342] + 7307872110 * arrold[344] + 3796297200 * arrold[346] +
       1251677700 * arrold[348] + 254186856 * arrold[350] + 30260340 * arrold[352] + 1947792 * arrold[354] +
       58905 * arrold[356] + 630 * arrold[358] + arrold[360]) * coeff;
#endif


#if EXPANSION > 18
      coeff = 1 / 262144.0;
      arr[361] =
      (arrold[361] + 703 * arrold[363] + 73815 * arrold[365] + 2760681 * arrold[367] +
       48903492 * arrold[369] + 472733756 * arrold[371] + 2707475148 * arrold[373] +
       9669554100 * arrold[375] + 22239974430 * arrold[377] + 33578000610 * arrold[379]) * coeff;
      arr[362] =
      (36 * arrold[362] + 7104 * arrold[364] + 369852 * arrold[366] + 7970688 * arrold[368] +
       85795600 * arrold[370] + 506662016 * arrold[372] + 1709984304 * arrold[374] +
       3257112960 * arrold[376] + 3029594040 * arrold[378]) * coeff;
      arr[363] =
      (arrold[361] + 559 * arrold[363] + 45255 * arrold[365] + 1252713 * arrold[367] +
       15512772 * arrold[369] + 96160636 * arrold[371] + 304253964 * arrold[373] + 426395700 * arrold[375] -
       31635810 * arrold[377] - 811985790 * arrold[379]) * coeff;
      arr[364] =
      (32 * arrold[362] + 4864 * arrold[364] + 186592 * arrold[366] + 2776576 * arrold[368] +
       18550400 * arrold[370] + 54774272 * arrold[372] + 41080704 * arrold[374] - 117373440 * arrold[376] -
       218349120 * arrold[378]) * coeff;
      arr[365] =
      (arrold[361] + 431 * arrold[363] + 25671 * arrold[365] + 486761 * arrold[367] + 3640516 * arrold[369] +
       10086780 * arrold[371] - 916980 * arrold[373] - 43098060 * arrold[375] - 31635810 * arrold[377] +
       61410690 * arrold[379]) * coeff;
      arr[366] =
      (28 * arrold[362] + 3136 * arrold[364] + 82180 * arrold[366] + 725120 * arrold[368] +
       1936880 * arrold[370] - 2186368 * arrold[372] - 12212016 * arrold[374] + 1726080 * arrold[376] +
       27293640 * arrold[378]) * coeff;
      arr[367] =
      (arrold[361] + 319 * arrold[363] + 13015 * arrold[365] + 145385 * arrold[367] + 398660 * arrold[369] -
       902596 * arrold[371] - 3160884 * arrold[373] + 3506100 * arrold[375] + 8064030 * arrold[377] -
       8064030 * arrold[379]) * coeff;
      arr[368] =
      (24 * arrold[362] + 1856 * arrold[364] + 28840 * arrold[366] + 90240 * arrold[368] -
       292640 * arrold[370] - 805504 * arrold[372] + 1812384 * arrold[374] + 1726080 * arrold[376] -
       4962480 * arrold[378]) * coeff;
      arr[369] =
      (arrold[361] + 223 * arrold[363] + 5495 * arrold[365] + 22505 * arrold[367] - 85180 * arrold[369] -
       215876 * arrold[371] + 747852 * arrold[373] + 165300 * arrold[375] - 2181090 * arrold[377] +
       1540770 * arrold[379]) * coeff;
      arr[370] =
      (20 * arrold[362] + 960 * arrold[364] + 5964 * arrold[366] - 22656 * arrold[368] - 64816 * arrold[370] +
       285824 * arrold[372] - 87696 * arrold[374] - 835200 * arrold[376] + 1200600 * arrold[378]) * coeff;
      arr[371] =
      (arrold[361] + 143 * arrold[363] + 1575 * arrold[365] - 5271 * arrold[367] - 22332 * arrold[369] +
       106236 * arrold[371] - 73332 * arrold[373] - 305100 * arrold[375] + 689310 * arrold[377] -
       391230 * arrold[379]) * coeff;
      arr[372] =
      (16 * arrold[362] + 384 * arrold[364] - 912 * arrold[366] - 8448 * arrold[368] + 38720 * arrold[370] -
       35584 * arrold[372] - 115776 * arrold[374] + 357120 * arrold[376] - 364320 * arrold[378]) * coeff;
      arr[373] =
      (arrold[361] + 79 * arrold[363] - 25 * arrold[365] - 3223 * arrold[367] + 13508 * arrold[369] -
       12804 * arrold[371] - 50036 * arrold[373] + 181300 * arrold[375] - 252770 * arrold[377] +
       123970 * arrold[379]) * coeff;
      arr[374] =
      (12 * arrold[362] + 64 * arrold[364] - 1132 * arrold[366] + 4224 * arrold[368] - 2640 * arrold[370] -
       25728 * arrold[372] + 94224 * arrold[374] - 158080 * arrold[376] + 131560 * arrold[378]) * coeff;
      arr[375] =
      (arrold[361] + 31 * arrold[363] - 329 * arrold[365] + 1001 * arrold[367] + 836 * arrold[369] -
       14916 * arrold[371] + 50764 * arrold[373] - 94796 * arrold[375] + 103454 * arrold[377] -
       46046 * arrold[379]) * coeff;
      arr[376] =
      (8 * arrold[362] - 64 * arrold[364] + 56 * arrold[366] + 1408 * arrold[368] - 8800 * arrold[370] +
       27776 * arrold[372] - 55328 * arrold[374] + 71552 * arrold[376] - 52624 * arrold[378]) * coeff;
      arr[377] =
      (arrold[361] - arrold[363] - 105 * arrold[365] + 1001 * arrold[367] - 4796 * arrold[369] +
       14652 * arrold[371] - 30772 * arrold[373] + 44980 * arrold[375] - 42978 * arrold[377] +
       18018 * arrold[379]) * coeff;
      arr[378] =
      (4 * arrold[362] - 64 * arrold[364] + 476 * arrold[366] - 2176 * arrold[368] + 6800 * arrold[370] -
       15232 * arrold[372] + 24752 * arrold[374] - 28288 * arrold[376] + 19448 * arrold[378]) * coeff;
      arr[379] =
      (arrold[361] - 17 * arrold[363] + 135 * arrold[365] - 663 * arrold[367] + 2244 * arrold[369] -
       5508 * arrold[371] + 9996 * arrold[373] - 13260 * arrold[375] + 11934 * arrold[377] -
       4862 * arrold[379]) * coeff;
      arr[380] =
      (-92378 * arrold[381] + 75582 * arrold[383] - 50388 * arrold[385] + 27132 * arrold[387] -
       11628 * arrold[389] + 3876 * arrold[391] - 969 * arrold[393] + 171 * arrold[395] - 19 * arrold[397] +
       arrold[399]) * coeff;
      arr[381] =
      (-48620 * arrold[380] + 87516 * arrold[382] - 63648 * arrold[384] + 37128 * arrold[386] -
       17136 * arrold[388] + 6120 * arrold[390] - 1632 * arrold[392] + 306 * arrold[394] - 36 * arrold[396] +
       2 * arrold[398]) * coeff;
      arr[382] =
      (102102 * arrold[381] - 80002 * arrold[383] + 48620 * arrold[385] - 22372 * arrold[387] +
       7412 * arrold[389] - 1564 * arrold[391] + 119 * arrold[393] + 35 * arrold[395] - 11 * arrold[397] +
       arrold[399]) * coeff;
      arr[383] =
      (60060 * arrold[380] - 103532 * arrold[382] + 65312 * arrold[384] - 28392 * arrold[386] +
       6832 * arrold[388] + 440 * arrold[390] - 1056 * arrold[392] + 406 * arrold[394] - 76 * arrold[396] +
       6 * arrold[398]) * coeff;
      arr[384] =
      (-138138 * arrold[381] + 93886 * arrold[383] - 38740 * arrold[385] + 3836 * arrold[387] +
       6292 * arrold[389] - 4444 * arrold[391] + 1463 * arrold[393] - 245 * arrold[395] + 13 * arrold[397] +
       arrold[399]) * coeff;
      arr[385] =
      (-92092 * arrold[380] + 144716 * arrold[382] - 61984 * arrold[384] - 728 * arrold[386] +
       19152 * arrold[388] - 12408 * arrold[390] + 3872 * arrold[392] - 518 * arrold[394] - 20 * arrold[396] +
       10 * arrold[398]) * coeff;
      arr[386] =
      (230230 * arrold[381] - 116610 * arrold[383] - 1300 * arrold[385] + 44188 * arrold[387] -
       29964 * arrold[389] + 8932 * arrold[391] - 649 * arrold[393] - 285 * arrold[395] + 53 * arrold[397] +
       arrold[399]) * coeff;
      arr[387] =
      (177100 * arrold[380] - 237820 * arrold[382] + 21920 * arrold[384] + 88376 * arrold[386] -
       68496 * arrold[388] + 19800 * arrold[390] + 352 * arrold[392] - 1442 * arrold[394] +
       196 * arrold[396] + 14 * arrold[398]) * coeff;
      arr[388] =
      (-478170 * arrold[381] + 126270 * arrold[383] + 153900 * arrold[385] - 154116 * arrold[387] +
       45716 * arrold[389] + 5412 * arrold[391] - 5577 * arrold[393] + 555 * arrold[395] + 109 * arrold[397] +
       arrold[399]) * coeff;
      arr[389] =
      (-434700 * arrold[380] + 451260 * arrold[382] + 205920 * arrold[384] - 343224 * arrold[386] +
       116368 * arrold[388] + 21800 * arrold[390] - 19296 * arrold[392] + 1218 * arrold[394] +
       636 * arrold[396] + 18 * arrold[398]) * coeff;
      arr[390] =
      (1260630 * arrold[381] + 60030 * arrold[383] - 736020 * arrold[385] + 328860 * arrold[387] +
       63220 * arrold[389] - 64284 * arrold[391] + 1911 * arrold[393] + 3171 * arrold[395] +
       181 * arrold[397] + arrold[399]) * coeff;
      arr[391] =
      (1400700 * arrold[380] - 920460 * arrold[382] - 1405920 * arrold[384] + 989016 * arrold[386] +
       133168 * arrold[388] - 214280 * arrold[390] + 1760 * arrold[392] + 14630 * arrold[394] +
       1364 * arrold[396] + 22 * arrold[398]) * coeff;
      arr[392] =
      (-4342170 * arrold[381] - 1860930 * arrold[383] + 2966700 * arrold[385] + 75516 * arrold[387] -
       722796 * arrold[389] + 6820 * arrold[391] + 65975 * arrold[393] + 8715 * arrold[395] +
       269 * arrold[397] + arrold[399]) * coeff;
      arr[393] =
      (-6203100 * arrold[380] + 1240620 * arrold[382] + 8198880 * arrold[384] - 1272984 * arrold[386] -
       2430896 * arrold[388] + 112840 * arrold[390] + 301600 * arrold[392] + 50570 * arrold[394] +
       2444 * arrold[396] + 26 * arrold[398]) * coeff;
      arr[394] =
      (20470230 * arrold[381] + 17988990 * arrold[383] - 9978900 * arrold[385] - 7778148 * arrold[387] +
       1147124 * arrold[389] + 1425380 * arrold[391] + 278135 * arrold[393] + 18595 * arrold[395] +
       373 * arrold[397] + arrold[399]) * coeff;
      arr[395] =
      (40940460 * arrold[380] + 13646820 * arrold[382] - 51350880 * arrold[384] - 20907144 * arrold[386] +
       9047536 * arrold[388] + 7002776 * arrold[390] + 1490016 * arrold[392] + 126446 * arrold[394] +
       3940 * arrold[396] + 30 * arrold[398]) * coeff;
      arr[396] =
      (-143291610 * arrold[381] - 200360130 * arrold[383] - 22924500 * arrold[385] + 62904828 * arrold[387] +
       35639956 * arrold[389] + 7907108 * arrold[391] + 799799 * arrold[393] + 34475 * arrold[395] +
       493 * arrold[397] + arrold[399]) * coeff;
      arr[397] =
      (-477638700 * arrold[380] - 450345060 * arrold[382] + 286097760 * arrold[384] +
       408239496 * arrold[386] + 186574864 * arrold[388] + 41970280 * arrold[390] + 4829088 * arrold[392] +
       266322 * arrold[394] + 5916 * arrold[396] + 34 * arrold[398]) * coeff;
      arr[398] =
      (1767263190 * arrold[381] + 3511574910 * arrold[383] + 2544619500 * arrold[385] +
       997490844 * arrold[387] + 223926516 * arrold[389] + 28312548 * arrold[391] + 1888887 * arrold[393] +
       58275 * arrold[395] + 629 * arrold[397] + arrold[399]) * coeff;
      arr[399] =
      (17672631900 * arrold[380] + 28781143380 * arrold[382] + 15471286560 * arrold[384] +
       5414950296 * arrold[386] + 1203322288 * arrold[388] + 163011640 * arrold[390] +
       12620256 * arrold[392] + 501942 * arrold[394] + 8436 * arrold[396] + 38 * arrold[398]) * coeff;
#endif


#if EXPANSION > 19
      coeff = 1 / 524288.0;
      arr[400] =
      (40 * arrold[401] + 9880 * arrold[403] + 658008 * arrold[405] + 18643560 * arrold[407] +
       273438880 * arrold[409] + 2311801440 * arrold[411] + 12033222880 * arrold[413] +
       40225345056 * arrold[415] + 88732378800 * arrold[417] + 131282408400 * arrold[419]) * coeff;
      arr[401] =
      (arrold[400] + 702 * arrold[402] + 73112 * arrold[404] + 2686866 * arrold[406] +
       46142811 * arrold[408] + 423830264 * arrold[410] + 2234741392 * arrold[412] +
       6962078952 * arrold[414] + 12570420330 * arrold[416] + 11338026180 * arrold[418]) * coeff;
      arr[402] =
      (36 * arrold[401] + 7068 * arrold[403] + 362748 * arrold[405] + 7600836 * arrold[407] +
       77824912 * arrold[409] + 420866416 * arrold[411] + 1203322288 * arrold[413] +
       1547128656 * arrold[415] - 227518920 * arrold[417] - 3029594040 * arrold[419]) * coeff;
      arr[403] =
      (arrold[400] + 558 * arrold[402] + 44696 * arrold[404] + 1207458 * arrold[406] +
       14260059 * arrold[408] + 80647864 * arrold[410] + 208093328 * arrold[412] + 122141736 * arrold[414] -
       458031510 * arrold[416] - 780349980 * arrold[418]) * coeff;
      arr[404] =
      (32 * arrold[401] + 4832 * arrold[403] + 181728 * arrold[405] + 2589984 * arrold[407] +
       15773824 * arrold[409] + 36223872 * arrold[411] - 13693568 * arrold[413] - 158454144 * arrold[415] -
       100975680 * arrold[417] + 218349120 * arrold[419]) * coeff;
      arr[405] =
      (arrold[400] + 430 * arrold[402] + 25240 * arrold[404] + 461090 * arrold[406] + 3153755 * arrold[408] +
       6446264 * arrold[410] - 11003760 * arrold[412] - 42181080 * arrold[414] + 11462250 * arrold[416] +
       93046500 * arrold[418]) * coeff;
      arr[406] =
      (28 * arrold[401] + 3108 * arrold[403] + 79044 * arrold[405] + 642940 * arrold[407] +
       1211760 * arrold[409] - 4123248 * arrold[411] - 10025648 * arrold[413] + 13938096 * arrold[415] +
       25567560 * arrold[417] - 27293640 * arrold[419]) * coeff;
      arr[407] =
      (arrold[400] + 318 * arrold[402] + 12696 * arrold[404] + 132370 * arrold[406] + 253275 * arrold[408] -
       1301256 * arrold[410] - 2258288 * arrold[412] + 6666984 * arrold[414] + 4557930 * arrold[416] -
       16128060 * arrold[418]) * coeff;
      arr[408] =
      (24 * arrold[401] + 1832 * arrold[403] + 26984 * arrold[405] + 61400 * arrold[407] -
       382880 * arrold[409] - 512864 * arrold[411] + 2617888 * arrold[413] - 86304 * arrold[415] -
       6688560 * arrold[417] + 4962480 * arrold[419]) * coeff;
      arr[409] =
      (arrold[400] + 222 * arrold[402] + 5272 * arrold[404] + 17010 * arrold[406] - 107685 * arrold[408] -
       130696 * arrold[410] + 963728 * arrold[412] - 582552 * arrold[414] - 2346390 * arrold[416] +
       3721860 * arrold[418]) * coeff;
      arr[410] =
      (20 * arrold[401] + 940 * arrold[403] + 5004 * arrold[405] - 28620 * arrold[407] - 42160 * arrold[409] +
       350640 * arrold[411] - 373520 * arrold[413] - 747504 * arrold[415] + 2035800 * arrold[417] -
       1200600 * arrold[419]) * coeff;
      arr[411] =
      (arrold[400] + 142 * arrold[402] + 1432 * arrold[404] - 6846 * arrold[406] - 17061 * arrold[408] +
       128568 * arrold[410] - 179568 * arrold[412] - 231768 * arrold[414] + 994410 * arrold[416] -
       1080540 * arrold[418]) * coeff;
      arr[412] =
      (16 * arrold[401] + 368 * arrold[403] - 1296 * arrold[405] - 7536 * arrold[407] + 47168 * arrold[409] -
       74304 * arrold[411] - 80192 * arrold[413] + 472896 * arrold[415] - 721440 * arrold[417] +
       364320 * arrold[419]) * coeff;
      arr[413] =
      (arrold[400] + 78 * arrold[402] - 104 * arrold[404] - 3198 * arrold[406] + 16731 * arrold[408] -
       26312 * arrold[410] - 37232 * arrold[412] + 231336 * arrold[414] - 434070 * arrold[416] +
       376740 * arrold[418]) * coeff;
      arr[414] =
      (12 * arrold[401] + 52 * arrold[403] - 1196 * arrold[405] + 5356 * arrold[407] - 6864 * arrold[409] -
       23088 * arrold[411] + 119952 * arrold[413] - 252304 * arrold[415] + 289640 * arrold[417] -
       131560 * arrold[419]) * coeff;
      arr[415] =
      (arrold[400] + 30 * arrold[402] - 360 * arrold[404] + 1330 * arrold[406] - 165 * arrold[408] -
       15752 * arrold[410] + 65680 * arrold[412] - 145560 * arrold[414] + 198250 * arrold[416] -
       149500 * arrold[418]) * coeff;
      arr[416] =
      (8 * arrold[401] - 72 * arrold[403] + 120 * arrold[405] + 1352 * arrold[407] - 10208 * arrold[409] +
       36576 * arrold[411] - 83104 * arrold[413] + 126880 * arrold[415] - 124176 * arrold[417] +
       52624 * arrold[419]) * coeff;
      arr[417] =
      (arrold[400] - 2 * arrold[402] - 104 * arrold[404] + 1106 * arrold[406] - 5797 * arrold[408] +
       19448 * arrold[410] - 45424 * arrold[412] + 75752 * arrold[414] - 87958 * arrold[416] +
       60996 * arrold[418]) * coeff;
      arr[418] =
      (4 * arrold[401] - 68 * arrold[403] + 540 * arrold[405] - 2652 * arrold[407] + 8976 * arrold[409] -
       22032 * arrold[411] + 39984 * arrold[413] - 53040 * arrold[415] + 47736 * arrold[417] -
       19448 * arrold[419]) * coeff;
      arr[419] =
      (arrold[400] - 18 * arrold[402] + 152 * arrold[404] - 798 * arrold[406] + 2907 * arrold[408] -
       7752 * arrold[410] + 15504 * arrold[412] - 23256 * arrold[414] + 25194 * arrold[416] -
       16796 * arrold[418]) * coeff;
      arr[420] =
      (92378 * arrold[420] - 167960 * arrold[422] + 125970 * arrold[424] - 77520 * arrold[426] +
       38760 * arrold[428] - 15504 * arrold[430] + 4845 * arrold[432] - 1140 * arrold[434] +
       190 * arrold[436] - 20 * arrold[438] + arrold[440]) * coeff;
      arr[421] =
      (-184756 * arrold[421] + 151164 * arrold[423] - 100776 * arrold[425] + 54264 * arrold[427] -
       23256 * arrold[429] + 7752 * arrold[431] - 1938 * arrold[433] + 342 * arrold[435] - 38 * arrold[437] +
       2 * arrold[439]) * coeff;
      arr[422] =
      (-102102 * arrold[420] + 182104 * arrold[422] - 128622 * arrold[424] + 70992 * arrold[426] -
       29784 * arrold[428] + 8976 * arrold[430] - 1683 * arrold[432] + 84 * arrold[434] + 46 * arrold[436] -
       12 * arrold[438] + arrold[440]) * coeff;
      arr[423] =
      (223652 * arrold[421] - 168844 * arrold[423] + 93704 * arrold[425] - 35224 * arrold[427] +
       6392 * arrold[429] + 1496 * arrold[431] - 1462 * arrold[433] + 482 * arrold[435] - 82 * arrold[437] +
       6 * arrold[439]) * coeff;
      arr[424] =
      (138138 * arrold[420] - 232024 * arrold[422] + 132626 * arrold[424] - 42576 * arrold[426] -
       2456 * arrold[428] + 10736 * arrold[430] - 5907 * arrold[432] + 1708 * arrold[434] -
       258 * arrold[436] + 12 * arrold[438] + arrold[440]) * coeff;
      arr[425] =
      (-328900 * arrold[421] + 206700 * arrold[423] - 61256 * arrold[425] - 19880 * arrold[427] +
       31560 * arrold[429] - 16280 * arrold[431] + 4390 * arrold[433] - 498 * arrold[435] - 30 * arrold[437] +
       10 * arrold[439]) * coeff;
      arr[426] =
      (-230230 * arrold[420] + 346840 * arrold[422] - 115310 * arrold[424] - 45488 * arrold[426] +
       74152 * arrold[428] - 38896 * arrold[430] + 9581 * arrold[432] - 364 * arrold[434] -
       338 * arrold[436] + 52 * arrold[438] + arrold[440]) * coeff;
      arr[427] =
      (592020 * arrold[421] - 259740 * arrold[423] - 66456 * arrold[425] + 156872 * arrold[427] -
       88296 * arrold[429] + 19448 * arrold[431] + 1794 * arrold[433] - 1638 * arrold[435] +
       182 * arrold[437] + 14 * arrold[439]) * coeff;
      arr[428] =
      (478170 * arrold[420] - 604440 * arrold[422] - 27630 * arrold[424] + 308016 * arrold[426] -
       199832 * arrold[428] + 40304 * arrold[430] + 10989 * arrold[432] - 6132 * arrold[434] +
       446 * arrold[436] + 108 * arrold[438] + arrold[440]) * coeff;
      arr[429] =
      (-1320660 * arrold[421] + 245340 * arrold[423] + 549144 * arrold[425] - 459592 * arrold[427] +
       94568 * arrold[429] + 41096 * arrold[431] - 20514 * arrold[433] + 582 * arrold[435] +
       618 * arrold[437] + 18 * arrold[439]) * coeff;
      arr[430] =
      (-1260630 * arrold[420] + 1200600 * arrold[422] + 796050 * arrold[424] - 1064880 * arrold[426] +
       265640 * arrold[428] + 127504 * arrold[430] - 66195 * arrold[432] - 1260 * arrold[434] +
       2990 * arrold[436] + 180 * arrold[438] + arrold[440]) * coeff;
      arr[431] =
      (3721860 * arrold[421] + 485460 * arrold[423] - 2394936 * arrold[425] + 855848 * arrold[427] +
       347448 * arrold[429] - 216040 * arrold[431] - 12870 * arrold[433] + 13266 * arrold[435] +
       1342 * arrold[437] + 22 * arrold[439]) * coeff;
      arr[432] =
      (4342170 * arrold[420] - 2481240 * arrold[422] - 4827630 * arrold[424] + 2891184 * arrold[426] +
       798312 * arrold[428] - 729616 * arrold[430] - 59155 * arrold[432] + 57260 * arrold[434] +
       8446 * arrold[436] + 268 * arrold[438] + arrold[440]) * coeff;
      arr[433] =
      (-13646820 * arrold[421] - 6958260 * arrold[423] + 9471864 * arrold[425] + 1157912 * arrold[427] -
       2543736 * arrold[429] - 188760 * arrold[431] + 251030 * arrold[433] + 48126 * arrold[435] +
       2418 * arrold[437] + 26 * arrold[439]) * coeff;
      arr[434] =
      (-20470230 * arrold[420] + 2481240 * arrold[422] + 27967890 * arrold[424] - 2200752 * arrold[426] -
       8925272 * arrold[428] - 278256 * arrold[430] + 1147245 * arrold[432] + 259540 * arrold[434] +
       18222 * arrold[436] + 372 * arrold[438] + arrold[440]) * coeff;
      arr[435] =
      (68234100 * arrold[421] + 64997700 * arrold[423] - 30443736 * arrold[425] - 29954680 * arrold[427] +
       2044760 * arrold[429] + 5512760 * arrold[431] + 1363570 * arrold[433] + 122506 * arrold[435] +
       3910 * arrold[437] + 30 * arrold[439]) * coeff;
      arr[436] =
      (143291610 * arrold[420] + 57068520 * arrold[422] - 177435630 * arrold[424] - 85829328 * arrold[426] +
       27264872 * arrold[428] + 27732848 * arrold[430] + 7107309 * arrold[432] + 765324 * arrold[434] +
       33982 * arrold[436] + 492 * arrold[438] + arrold[440]) * coeff;
      arr[437] =
      (-504932340 * arrold[421] - 736442820 * arrold[423] - 122141736 * arrold[425] +
       221664632 * arrold[427] + 144604584 * arrold[429] + 37141192 * arrold[431] + 4562766 * arrold[433] +
       260406 * arrold[435] + 5882 * arrold[437] + 34 * arrold[439]) * coeff;
      arr[438] =
      (-1767263190 * arrold[420] - 1744311720 * arrold[422] + 966955410 * arrold[424] +
       1547128656 * arrold[426] + 773564328 * arrold[428] + 195613968 * arrold[430] + 26423661 * arrold[432] +
       1830612 * arrold[434] + 57646 * arrold[436] + 628 * arrold[438] + arrold[440]) * coeff;
      arr[439] =
      (6564120420 * arrold[421] + 13309856820 * arrold[423] + 10056336264 * arrold[425] +
       4211628008 * arrold[427] + 1040310648 * arrold[429] + 150391384 * arrold[431] +
       12118314 * arrold[433] + 493506 * arrold[435] + 8398 * arrold[437] + 38 * arrold[439]) * coeff;
      arr[440] =
      (68923264410 * arrold[420] + 113380261800 * arrold[422] + 62852101650 * arrold[424] +
       23206929840 * arrold[426] + 5586853480 * arrold[428] + 847660528 * arrold[430] +
       76904685 * arrold[432] + 3838380 * arrold[434] + 91390 * arrold[436] + 780 * arrold[438] +
       arrold[440]) * coeff;
#endif
  }
}
