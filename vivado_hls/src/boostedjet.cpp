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

void boostedjet(ap_uint<10> jet_seed,             // input
			  region_t regions[NR_CNTR_REG],
			  ap_uint<10> et_3by3[NR_CNTR_REG], // input 14x18
			  ap_uint<10> et_jet [NR_SCNTR_REG], // *output* 7x9
			  bitset<3> rEta_jet [NR_SCNTR_REG],
			  bitset<3> rPhi_jet [NR_SCNTR_REG])
 
{

	bool jet_veto[NR_CNTR_REG];
	ap_uint<10> sr_et[NR_SCNTR_REG];
	bitset<3> sr_eta[NR_SCNTR_REG];
	bitset<3> sr_phi[NR_SCNTR_REG];
	bool activeRegion[9];

#pragma HLS INTERFACE ap_none port=jet_seed

#pragma HLS PIPELINE II=3 // target clk freq = 120 MHz

#pragma HLS ARRAY_RESHAPE  variable=regions   complete  dim=1
#pragma HLS ARRAY_RESHAPE  variable=et_3by3   complete  dim=1
#pragma HLS ARRAY_RESHAPE  variable=et_jet    complete  dim=1
#pragma HLS ARRAY_RESHAPE  variable=rEta_jet  complete  dim=1
#pragma HLS ARRAY_RESHAPE  variable=rPhi_jet  complete  dim=1
#pragma HLS ARRAY_RESHAPE  variable=jet_veto  complete  dim=1
#pragma HLS ARRAY_RESHAPE  variable=activeRegion complete  dim=1
#pragma HLS ARRAY_PARTITION  variable=sr_et     complete  dim=1
#pragma HLS ARRAY_PARTITION  variable=sr_eta    complete  dim=1
#pragma HLS ARRAY_PARTITION  variable=sr_phi    complete  dim=1
#pragma HLS inline region

	for (int idx = 0; idx < NR_SCNTR_REG; idx++)
	{
#pragma HLS UNROLL
		sr_et[idx] = 0;
		sr_eta[idx] = 0;
		sr_phi[idx] = 0;
	}

	loop_rgn_et: for (int idx = 0; idx < NR_CNTR_REG; idx++)
	{
#pragma HLS UNROLL

		ap_uint<8> sidx = central_super_region_idx[idx];

		if(!(idx % 14 == 0) && !((idx + 1) % 14 == 0))
		{
			ap_uint<10> et_C, et_N, et_S, et_E, et_W, et_NE, et_NW, et_SE, et_SW;
			bool tauveto_C, tauveto_N, tauveto_S, tauveto_E, tauveto_W, tauveto_NE, tauveto_NW, tauveto_SE, tauveto_SW;
			ap_int<10> neigh_N = region_cntr_neighbors[idx].nb_N;
			ap_int<10> neigh_S = region_cntr_neighbors[idx].nb_S;
			ap_int<10> neigh_E = region_cntr_neighbors[idx].nb_E;
			ap_int<10> neigh_W = region_cntr_neighbors[idx].nb_W;
			ap_int<10> neigh_NE = region_cntr_neighbors[idx].nb_NE;
			ap_int<10> neigh_NW = region_cntr_neighbors[idx].nb_NW;
			ap_int<10> neigh_SE = region_cntr_neighbors[idx].nb_SE;
			ap_int<10> neigh_SW = region_cntr_neighbors[idx].nb_SW;

			et_C = regions[idx].et; 
			tauveto_C = regions[idx].tau_veto;

			if (neigh_N != -1) {
				et_N = regions[neigh_N].et; 
				tauveto_N = regions[neigh_N].tau_veto;
			}
			else {
				et_N = 0;
				tauveto_N = false;
			}

			if (neigh_S != -1) {
				et_S = regions[neigh_S].et;
				tauveto_S = regions[neigh_S].tau_veto;
			}
			else {
				et_S = 0;
				tauveto_S = false;
			}

			if (neigh_E != -1) {
				et_E  = regions[neigh_E].et; 
				tauveto_E  = regions[neigh_E].tau_veto;
			}
			else {
				et_E  = 0;
				tauveto_E  = false;
			}

			if (neigh_W != -1) {
				et_W  = regions[neigh_W].et;
				tauveto_W = regions[neigh_W].tau_veto;
			}
			else {
				et_W = 0;
				tauveto_W = false;
			}

			if (neigh_NE != -1) {
				et_NE = regions[neigh_NE].et; 
				tauveto_NE = regions[neigh_NE].tau_veto;
			}
			else {
				et_NE = 0;
				tauveto_NE = false;
			}

			if (neigh_NW != -1) {
				et_NW = regions[neigh_NW].et; 
				tauveto_NW = regions[neigh_NW].tau_veto;
			}
			else {
				et_NW = 0;
				tauveto_NW = false;
			}

			if (neigh_SE != -1) {
				et_SE = regions[neigh_SE].et; 
				tauveto_SE = regions[neigh_SE].tau_veto;
			}
			else {
				et_SE = 0;
				tauveto_SE = false;
			}

			if (neigh_SW != -1) {
				et_SW = regions[neigh_SW].et; 
				tauveto_SW = regions[neigh_SW].tau_veto;
			}
			else {
				et_SW = 0;
				tauveto_SW = false;
			}

			jet_veto[idx] = false;

			if (et_C < jet_seed)  jet_veto[idx] = true;
			if (et_C < et_N)      jet_veto[idx] = true;
			if (et_C < et_E)      jet_veto[idx] = true;
			if (et_C < et_W)      jet_veto[idx] = true;
			if (et_C < et_S)      jet_veto[idx] = true;
			if (et_C < et_NW)     jet_veto[idx] = true;
			if (et_C < et_NE)     jet_veto[idx] = true;
			if (et_C < et_SE)     jet_veto[idx] = true;
			if (et_C < et_SW)     jet_veto[idx] = true;

			// assign et_jet and pattern
			if (jet_veto[idx] == false && et_3by3[idx] > sr_et[sidx]){
				sr_et[sidx] = et_3by3[idx];
				et_jet[sidx] = et_3by3[idx];
				if (!tauveto_C && et_C > 30 && et_C > (et_jet[sidx] >> 4)) activeRegion[4] = true;
				else activeRegion[4] = false;
				if (!tauveto_N && et_N > 30 && et_N > (et_jet[sidx] >> 4)) activeRegion[3] = true;
				else activeRegion[3] = false;
				if (!tauveto_E && et_E > 30 && et_E > (et_jet[sidx] >> 4)) activeRegion[1] = true;
				else activeRegion[1] = false;
				if (!tauveto_W && et_W > 30 && et_W > (et_jet[sidx] >> 4)) activeRegion[7] = true;
				else activeRegion[7] = false;
				if (!tauveto_S && et_S > 30 && et_S > (et_jet[sidx] >> 4)) activeRegion[5] = true;
				else activeRegion[5] = false;
				if (!tauveto_NW && et_NW > 30 && et_NW > (et_jet[sidx] >> 4)) activeRegion[6] = true;
				else activeRegion[6] = false;
				if (!tauveto_NE && et_NE > 30 && et_NE > (et_jet[sidx] >> 4)) activeRegion[0] = true;
				else activeRegion[0] = false;
				if (!tauveto_SW && et_SW > 30 && et_SW > (et_jet[sidx] >> 4)) activeRegion[8] = true;
				else activeRegion[8] = false;
				if (!tauveto_SE && et_SE > 30 && et_SE > (et_jet[sidx] >> 4)) activeRegion[2] = true;
				else activeRegion[2] = false;
				rEta_jet[sidx] = etapattern(activeRegion);
				rPhi_jet[sidx] = phipattern(activeRegion);
				sr_eta[sidx] = rEta_jet[sidx];
				sr_phi[sidx] = rPhi_jet[sidx];
			}

			else {
				et_jet[sidx] = sr_et[sidx];
				rEta_jet[sidx] = sr_eta[sidx];
				rPhi_jet[sidx] = sr_phi[sidx];
			}

		}
	}
	return;
}
