#include "ap_int.h"
#include "UCTSummaryCard.hpp"
#include "region_neighbors.h"

void superregion_centre(ap_uint<10> et[NR_CNTR_REG], int sr_centre_index[NR_SCNTR_REG])
{
#pragma HLS PIPELINE II=3
#pragma HLS ARRAY_RESHAPE variable=et complete dim=1
#pragma HLS ARRAY_RESHAPE variable=sr_centre_index complete dim=1

	int idx, sidx, neigh_N, neigh_E, neigh_NE;
	ap_uint<10> tmp, et_N, et_E, et_NE;

	scntr_phi_loop: for (unsigned int phi = 0; phi < 9; phi++){
#pragma HLS UNROLL
		scntr_reg_loop: for (int reg = 3; reg < 10; reg++){
#pragma HLS UNROLL
			idx = 28*phi - 6 + 2*reg;
			sidx = 7*phi - 3 + reg;
			tmp = et[idx];
			sr_centre_index[sidx] = idx;
			neigh_N = region_cntr_neighbors[idx].nb_N;
			neigh_E = region_cntr_neighbors[idx].nb_E;
			neigh_NE = region_cntr_neighbors[idx].nb_NE;
			if(neigh_N != -1)
				et_N = et[neigh_N];
			else
				et_N = 0;
			if(neigh_E != -1)
				et_E = et[neigh_E];
			else
				et_E = 0;
			if(neigh_NE != -1)
				et_NE = et[neigh_NE];
			else
				et_NE = 0;
			if(et_N > tmp) { sr_centre_index[sidx] = neigh_N; tmp = et_N; }
			if(et_E > tmp) { sr_centre_index[sidx] = neigh_E; tmp = et_E; }
			if(et_NE > tmp) { sr_centre_index[sidx] = neigh_NE; tmp = et_NE; }
		}
	}		

	return;
}

void superregion(ap_uint<10> et[NR_CALO_REG], int sr_index[NR_SUPER_REG])
{
#pragma HLS PIPELINE II=3
#pragma HLS ARRAY_RESHAPE variable=et complete dim=1
#pragma HLS ARRAY_RESHAPE variable=sr_index complete dim=1

	int idx, sidx, neigh_N, neigh_E, neigh_NE;
	ap_uint<10> tmp, et_N, et_E, et_NE;

	sregion_phi_loop: for (unsigned int phi = 0; phi < 9; phi++){
#pragma HLS UNROLL
		sregion_reg_loop: for (int reg = 0; reg < 13; reg++){
#pragma HLS UNROLL
			idx = 52*phi + 2*reg;
			sidx = 13*phi + reg;
			tmp = et[idx];
			sr_index[sidx] = idx;
			neigh_N = rgn_nghbr[idx].nb_N;
			neigh_E = rgn_nghbr[idx].nb_E;
			neigh_NE = rgn_nghbr[idx].nb_NE;
			if(neigh_N != -1)
				et_N = et[neigh_N];
			else
				et_N = 0;
			if(neigh_E != -1)
				et_E = et[neigh_E];
			else
				et_E = 0;
			if(neigh_NE != -1)
				et_NE = et[neigh_NE];
			else
				et_NE = 0;
			if(et_N > tmp) { sr_index[sidx] = neigh_N; tmp = et_N; }
			if(et_E > tmp) { sr_index[sidx] = neigh_E; tmp = et_E; }
			if(et_NE > tmp) { sr_index[sidx] = neigh_NE; tmp = et_NE; }
		}
	}

	return;
} 
