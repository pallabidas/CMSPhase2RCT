#ifndef __SUPERREGION_H__
#define __SUPERREGION_H__

#include "UCTSummaryCard.hpp"

const ap_uint<8> central_super_region_idx[NR_CNTR_REG] = { 0, 1, 1, 1, 1, 2, 2, 2, 2, 3, 3, 3, 3, 0, 0, 1, 1, 1, 1, 2, 2, 2, 2, 3, 3, 3, 3, 0, 0, 1, 1, 1, 1, 2, 2, 2, 2, 3, 3, 3, 3, 0, 0, 1, 1, 1, 1, 2, 2, 2, 2, 3, 3, 3, 3, 0, 0, 4, 4, 4, 4, 5, 5, 5, 5, 6, 6, 6, 6, 0, 0, 4, 4, 4, 4, 5, 5, 5, 5, 6, 6, 6, 6, 0, 0, 4, 4, 4, 4, 5, 5, 5, 5, 6, 6, 6, 6, 0, 0, 4, 4, 4, 4, 5, 5, 5, 5, 6, 6, 6, 6, 0, 0, 7, 7, 7, 7, 8, 8, 8, 8, 9, 9, 9, 9, 0, 0, 7, 7, 7, 7, 8, 8, 8, 8, 9, 9, 9, 9, 0, 0, 7, 7, 7, 7, 8, 8, 8, 8, 9, 9, 9, 9, 0, 0, 7, 7, 7, 7, 8, 8, 8, 8, 9, 9, 9, 9, 0, 0, 10, 10, 10, 10, 11, 11, 11, 11, 12, 12, 12, 12, 0, 0, 10, 10, 10, 10, 11, 11, 11, 11, 12, 12, 12, 12, 0, 0, 10, 10, 10, 10, 11, 11, 11, 11, 12, 12, 12, 12, 0, 0, 10, 10, 10, 10, 11, 11, 11, 11, 12, 12, 12, 12, 0, 0, 13, 13, 13, 13, 14, 14, 14, 14, 15, 15, 15, 15, 0, 0, 13, 13, 13, 13, 14, 14, 14, 14, 15, 15, 15, 15, 0 };


const ap_uint<8> super_region_idx[NR_CALO_REG] = { 0, 0, 2, 2, 2, 2, 2, 16, 16, 16, 16, 17, 17, 17, 17, 18, 18, 18, 18, 8, 8, 8, 8, 8, 14, 14, 0, 0, 2, 2, 2, 2, 2, 16, 16, 16, 16, 17, 17, 17, 17, 18, 18, 18, 18, 8, 8, 8, 8, 8, 14, 14, 0, 0, 2, 2, 2, 2, 2, 16, 16, 16, 16, 17, 17, 17, 17, 18, 18, 18, 18, 8, 8, 8, 8, 8, 14, 14, 0, 0, 3, 3, 3, 3, 3, 16, 16, 16, 16, 17, 17, 17, 17, 18, 18, 18, 18, 9, 9, 9, 9, 9, 14, 14, 0, 0, 3, 3, 3, 3, 3, 19, 19, 19, 19, 20, 20, 20, 20, 21, 21, 21, 21, 9, 9, 9, 9, 9, 14, 14, 0, 0, 3, 3, 3, 3, 3, 19, 19, 19, 19, 20, 20, 20, 20, 21, 21, 21, 21, 9, 9, 9, 9, 9, 14, 14, 0, 0, 4, 4, 4, 4, 4, 19, 19, 19, 19, 20, 20, 20, 20, 21, 21, 21, 21, 10, 10, 10, 10, 10, 14, 14, 0, 0, 4, 4, 4, 4, 4, 19, 19, 19, 19, 20, 20, 20, 20, 21, 21, 21, 21, 10, 10, 10, 10, 10, 14, 14, 0, 0, 4, 4, 4, 4, 4, 22, 22, 22, 22, 23, 23, 23, 23, 24, 24, 24, 24, 10, 10, 10, 10, 10, 14, 14, 1, 1, 5, 5, 5, 5, 5, 22, 22, 22, 22, 23, 23, 23, 23, 24, 24, 24, 24, 11, 11, 11, 11, 11, 15, 15, 1, 1, 5, 5, 5, 5, 5, 22, 22, 22, 22, 23, 23, 23, 23, 24, 24, 24, 24, 11, 11, 11, 11, 11, 15, 15, 1, 1, 5, 5, 5, 5, 5, 22, 22, 22, 22, 23, 23, 23, 23, 24, 24, 24, 24, 11, 11, 11, 11, 11, 15, 15, 1, 1, 6, 6, 6, 6, 6, 25, 25, 25, 25, 26, 26, 26, 26, 27, 27, 27, 27, 12, 12, 12, 12, 12, 15, 15, 1, 1, 6, 6, 6, 6, 6, 25, 25, 25, 25, 26, 26, 26, 26, 27, 27, 27, 27, 12, 12, 12, 12, 12, 15, 15, 1, 1, 6, 6, 6, 6, 6, 25, 25, 25, 25, 26, 26, 26, 26, 27, 27, 27, 27, 12, 12, 12, 12, 12, 15, 15, 1, 1, 7, 7, 7, 7, 7, 25, 25, 25, 25, 26, 26, 26, 26, 27, 27, 27, 27, 13, 13, 13, 13, 13, 15, 15, 1, 1, 7, 7, 7, 7, 7, 28, 28, 28, 29, 29, 29, 30, 30, 30, 31, 31, 31, 13, 13, 13, 13, 13, 15, 15, 1, 1, 7, 7, 7, 7, 7, 28, 28, 28, 29, 29, 29, 30, 30, 30, 31, 31, 31, 13, 13, 13, 13, 13, 15, 15 };

#endif
