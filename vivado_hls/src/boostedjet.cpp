#include <cstdlib>
#include "ap_int.h"
#include "UCTSummaryCard.hpp"
#include "region_neighbors.h"
#include <bitset>
#include "superregion.h"

using std::bitset;

bitset<3> etapattern(bool activeRegion[9])
{
#pragma HLS PIPELINE II=1
#pragma HLS ARRAY_RESHAPE  variable=activeRegion complete  dim=1
#pragma HLS inline

	bitset<3> rEta;
	for(int iEta = 0; iEta < 3; iEta++){
#pragma HLS UNROLL
		bool activeStrip = false;
		for(int iPhi = 0; iPhi < 3; iPhi++){
#pragma HLS UNROLL
			if(activeRegion[iEta*3+iPhi]) activeStrip = true;
		}
		if(activeStrip) rEta |= (0x1 << iEta);
	}
	return rEta;
}

bitset<3> phipattern(bool activeRegion[9])
{
#pragma HLS PIPELINE II=1
#pragma HLS ARRAY_RESHAPE  variable=activeRegion complete  dim=1
#pragma HLS inline

	bitset<3> rPhi;
	for(int iPhi = 0; iPhi < 3; iPhi++){
#pragma HLS UNROLL
		bool activeStrip = false;
		for(int iEta = 0; iEta < 3; iEta++){
#pragma HLS UNROLL
			if(activeRegion[iEta*3+iPhi]) activeStrip = true;
		}
		if(activeStrip) rPhi |= (0x1 << iPhi);
	}
	return rPhi;
}

void boostedjet(ap_uint<10> jet_seed,
			  region_t regions[NR_CNTR_REG],
			  ap_uint<10> et_3by3[NR_CNTR_REG],
			  ap_uint<10> et_jet [NR_SCNTR_REG],
			  bitset<3> rEta_jet [NR_SCNTR_REG],
			  bitset<3> rPhi_jet [NR_SCNTR_REG],
			  ap_uint<9> rIdx_boostedjet[NR_SCNTR_REG])
 
{

	bool jet_veto[NR_CNTR_REG];
	ap_uint<10> sr_et[NR_SCNTR_REG];
	ap_uint<10> sr_et_jet[NR_SCNTR_REG];
	bitset<3> sr_eta[NR_SCNTR_REG];
	bitset<3> sr_phi[NR_SCNTR_REG];
	ap_uint<9> sr_idx[NR_SCNTR_REG];
	bool activeRegion[9];
	bitset<3> b1 = 0b010;
	bitset<3> b2 = 0b011;
	bitset<3> b3 = 0b110;

#pragma HLS INTERFACE ap_none port=jet_seed

#pragma HLS PIPELINE II=1

#pragma HLS ARRAY_RESHAPE  variable=regions   complete  dim=1
#pragma HLS ARRAY_RESHAPE  variable=et_3by3   complete  dim=1
#pragma HLS ARRAY_RESHAPE  variable=et_jet    complete  dim=1
#pragma HLS ARRAY_RESHAPE  variable=rEta_jet  complete  dim=1
#pragma HLS ARRAY_RESHAPE  variable=rPhi_jet  complete  dim=1
#pragma HLS ARRAY_RESHAPE  variable=rIdx_boostedjet      complete  dim=1
#pragma HLS ARRAY_RESHAPE  variable=jet_veto  complete  dim=1
#pragma HLS ARRAY_RESHAPE  variable=activeRegion complete  dim=1
#pragma HLS ARRAY_PARTITION variable=sr_et     complete  dim=1
#pragma HLS ARRAY_PARTITION variable=sr_et_jet complete  dim=1
#pragma HLS ARRAY_PARTITION variable=sr_eta    complete  dim=1
#pragma HLS ARRAY_PARTITION variable=sr_phi    complete  dim=1
#pragma HLS ARRAY_PARTITION variable=sr_idx    complete  dim=1
#pragma HLS inline

	for (int idx = 0; idx < NR_SCNTR_REG; idx++)
	{
#pragma HLS UNROLL
		sr_et[idx] = 0;
		sr_et_jet[idx] = 0;
		sr_eta[idx] = 0;
		sr_phi[idx] = 0;
		sr_idx[idx] = 0;
	}

	loop_rgn_et: for (int idx = 0; idx < NR_CNTR_REG; idx++)
	{
#pragma HLS UNROLL

		ap_uint<8> sidx = central_super_region_idx[idx];
		ap_uint<10> et_j;
		bitset<3> reta_jet;
		bitset<3> rphi_jet;

		if ((idx % 14 == 0) || ((idx + 1) % 14) == 0)
		{
			et_j = 0;
			reta_jet = 0;
			rphi_jet = 0;
		}
		else
		{
			ap_uint<10> et_C, et_N, et_S, et_E, et_W, et_NE, et_NW, et_SE, et_SW, et_f;
			ap_int<10> neigh_N = region_cntr_neighbors[idx].nb_N;
			ap_int<10> neigh_S = region_cntr_neighbors[idx].nb_S;
			ap_int<10> neigh_E = region_cntr_neighbors[idx].nb_E;
			ap_int<10> neigh_W = region_cntr_neighbors[idx].nb_W;
			ap_int<10> neigh_NE = region_cntr_neighbors[idx].nb_NE;
			ap_int<10> neigh_NW = region_cntr_neighbors[idx].nb_NW;
			ap_int<10> neigh_SE = region_cntr_neighbors[idx].nb_SE;
			ap_int<10> neigh_SW = region_cntr_neighbors[idx].nb_SW;

			et_C = regions[idx].et; 

			if (neigh_N != -1) {
				et_N = regions[neigh_N].et; 
			}
			else {
				et_N = 0;
			}

			if (neigh_S != -1) {
				et_S = regions[neigh_S].et;
			}
			else {
				et_S = 0;
			}

			if (neigh_E != -1) {
				et_E  = regions[neigh_E].et; 
			}
			else {
				et_E  = 0;
			}

			if (neigh_W != -1) {
				et_W  = regions[neigh_W].et;
			}
			else {
				et_W = 0;
			}

			if (neigh_NE != -1) {
				et_NE = regions[neigh_NE].et; 
			}
			else {
				et_NE = 0;
			}

			if (neigh_NW != -1) {
				et_NW = regions[neigh_NW].et; 
			}
			else {
				et_NW = 0;
			}

			if (neigh_SE != -1) {
				et_SE = regions[neigh_SE].et; 
			}
			else {
				et_SE = 0;
			}

			if (neigh_SW != -1) {
				et_SW = regions[neigh_SW].et; 
			}
			else {
				et_SW = 0;
			}

			if (et_C < jet_seed)  jet_veto[idx] = true;
			else if (et_C < et_N)      jet_veto[idx] = true;
			else if (et_C < et_E)      jet_veto[idx] = true;
			else if (et_C < et_W)      jet_veto[idx] = true;
			else if (et_C < et_S)      jet_veto[idx] = true;
			else if (et_C < et_NW)     jet_veto[idx] = true;
			else if (et_C < et_NE)     jet_veto[idx] = true;
			else if (et_C < et_SE)     jet_veto[idx] = true;
			else if (et_C < et_SW)     jet_veto[idx] = true;
			else jet_veto[idx] = false;

			// assign et_jet and pattern
			if (jet_veto[idx] == false)
			{
				et_j = et_C + et_N + et_E + et_W + et_S + et_NW + et_NE + et_SW + et_SE;
				et_f = et_j >> 4;
				if(et_C > 30 && et_C > et_f) activeRegion[4] = true;
				else activeRegion[4] = false;
				if(et_N > 30 && et_N > et_f) activeRegion[3] = true;
				else activeRegion[3] = false;
				if(et_E > 30 && et_E > et_f) activeRegion[1] = true;
				else activeRegion[1] = false;
				if(et_W > 30 && et_W > et_f) activeRegion[7] = true;
				else activeRegion[7] = false;
				if(et_S > 30 && et_S > et_f) activeRegion[5] = true;
				else activeRegion[5] = false;
				if(et_NW > 30 && et_NW > et_f) activeRegion[6] = true;
				else activeRegion[6] = false;
				if(et_NE > 30 && et_NE > et_f) activeRegion[0] = true;
				else activeRegion[0] = false;
				if(et_SW > 30 && et_SW > et_f) activeRegion[8] = true;
				else activeRegion[8] = false;
				if(et_SE > 30 && et_SE > et_f) activeRegion[2] = true;
				else activeRegion[2] = false;
				reta_jet = etapattern(activeRegion);
				rphi_jet = phipattern(activeRegion);
			}
			else
			{
				et_j = 0;
				et_f = 0;
				activeRegion[4] = false;
				activeRegion[3] = false;
				activeRegion[1] = false;
				activeRegion[7] = false;
				activeRegion[5] = false;
				activeRegion[6] = false;
				activeRegion[0] = false;
				activeRegion[8] = false;
				activeRegion[2] = false;
				reta_jet = 0;
				rphi_jet = 0;
			}
		}

		if (regions[idx].et > sr_et[sidx] && ((reta_jet == b1) || (reta_jet == b2) || (reta_jet == b3) || (rphi_jet == b1) || (rphi_jet == b2) || (rphi_jet == b3)))
		{
			sr_et[sidx] = regions[idx].et;
			sr_et_jet[sidx] = et_3by3[idx];
			sr_eta[sidx] = reta_jet;
			sr_phi[sidx] = rphi_jet;
			sr_idx[sidx] = idx;
		}
	}

	for (int idx = 0; idx < NR_SCNTR_REG; idx++)
	{
#pragma HLS UNROLL
		et_jet[idx] = sr_et_jet[idx];
		rEta_jet[idx] = sr_eta[idx];
		rPhi_jet[idx] = sr_phi[idx];
		rIdx_boostedjet[idx] = sr_idx[idx];
	}

	return;
}
