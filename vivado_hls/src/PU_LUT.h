#ifndef PU_LUT_H
#define PU_LUT_H

#include "UCTSummaryCard.hpp"

void pu_lut(ap_uint<5> pum_bin, ap_uint<10> et[NR_CALO_REG],
		ap_uint<10> pu_sub_et[NR_CALO_REG]);

void pu_lut_cntr(ap_uint<5> pum_bin, ap_uint<10> et[NR_CNTR_REG],
                ap_uint<10> pu_sub_et[NR_CNTR_REG]);

#endif
