#include <cstdlib>
#include "ap_int.h"
#include "UCTSummaryCard.hpp"
#include "region_neighbors.h"

ap_uint<3> etapattern(bool activeRegion[9])
{
#pragma HLS PIPELINE II=6
#pragma HLS ARRAY_RESHAPE  variable=activeRegion complete  dim=1

	ap_uint<3> rEta;
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

ap_uint<3> phipattern(bool activeRegion[9])
{
#pragma HLS PIPELINE II=6
#pragma HLS ARRAY_RESHAPE  variable=activeRegion complete  dim=1

	ap_uint<3> rPhi;
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
			  //ap_uint<10> et_rgn [NR_CNTR_REG], // input 26x18
			  ap_uint<10> et_3by3[NR_CNTR_REG], // input 26x18
			  ap_uint<10> et_jet [NR_CNTR_REG], // *output* 26x18
			  ap_uint<3> rEta_jet [NR_CNTR_REG],
			  ap_uint<3> rPhi_jet [NR_CNTR_REG])
 
{

	bool jet_veto[NR_CNTR_REG];
	bool activeRegion[9];

#pragma HLS INTERFACE ap_none port=jet_seed

#pragma HLS PIPELINE II=6 // target clk freq = 250 MHz

#pragma HLS ARRAY_RESHAPE  variable=regions   complete  dim=1
//#pragma HLS ARRAY_RESHAPE  variable=et_rgn    complete  dim=1
#pragma HLS ARRAY_RESHAPE  variable=et_3by3   complete  dim=1
#pragma HLS ARRAY_RESHAPE  variable=et_jet    complete  dim=1
#pragma HLS ARRAY_RESHAPE  variable=rEta_jet  complete  dim=1
#pragma HLS ARRAY_RESHAPE  variable=rPhi_jet  complete  dim=1
#pragma HLS ARRAY_RESHAPE  variable=jet_veto  complete  dim=1
#pragma HLS ARRAY_RESHAPE  variable=activeRegion complete  dim=1

	loop_rgn_et: for (int idx = 0; idx < NR_CNTR_REG; idx++)
	{
#pragma HLS UNROLL

		if ((idx % 14 == 0) || ((idx + 1) % 14) == 0)
		{
		        et_jet[idx] = 0;
		        rEta_jet[idx] = 0;
			rPhi_jet[idx] = 0;
		}
		else
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
			if (jet_veto[idx] == false){
				et_jet[idx] = et_3by3[idx];
				if (tauveto_C && et_C > (et_jet[idx] >> 4)) activeRegion[4] = true;
				else activeRegion[4] = false;
				if (tauveto_N && et_N > (et_jet[idx] >> 4)) activeRegion[3] = true;
				else activeRegion[3] = false;
				if (tauveto_E && et_E > (et_jet[idx] >> 4)) activeRegion[1] = true;
				else activeRegion[1] = false;
				if (tauveto_W && et_W > (et_jet[idx] >> 4)) activeRegion[7] = true;
				else activeRegion[7] = false;
				if (tauveto_S && et_S > (et_jet[idx] >> 4)) activeRegion[5] = true;
				else activeRegion[5] = false;
				if (tauveto_NW && et_NW > (et_jet[idx] >> 4)) activeRegion[6] = true;
				else activeRegion[6] = false;
				if (tauveto_NE && et_NE > (et_jet[idx] >> 4)) activeRegion[0] = true;
				else activeRegion[0] = false;
				if (tauveto_SW && et_SW > (et_jet[idx] >> 4)) activeRegion[8] = true;
				else activeRegion[8] = false;
				if (tauveto_SE && et_SE > (et_jet[idx] >> 4)) activeRegion[2] = true;
				else activeRegion[2] = false;
				rEta_jet[idx] = etapattern(activeRegion);
				rPhi_jet[idx] = phipattern(activeRegion);			
			}
			else {
				et_jet[idx] = 0;
				rEta_jet[idx] = 0;
				rPhi_jet[idx] = 0;
			}
		}

//		rEta_jet[idx] = 0;
//		rPhi_jet[idx] = 0;
//
//		if (regions[rgn_nghbr[idx].nb_C].tau_veto == true && et_C > (et_jet[idx] >> 4)) activeRegion[4] = true;
//		else activeRegion[4] = false;
//		if (regions[rgn_nghbr[idx].nb_N].tau_veto == true && et_N > (et_jet[idx] >> 4)) activeRegion[3] = true;
//		else activeRegion[3] = false;
//		if (regions[rgn_nghbr[idx].nb_E].tau_veto == true && et_E > (et_jet[idx] >> 4)) activeRegion[1] = true;
//		else activeRegion[1] = false;
//		if (regions[rgn_nghbr[idx].nb_W].tau_veto == true && et_W > (et_jet[idx] >> 4)) activeRegion[7] = true;
//		else activeRegion[7] = false;
//		if (regions[rgn_nghbr[idx].nb_S].tau_veto == true && et_S > (et_jet[idx] >> 4)) activeRegion[5] = true;
//		else activeRegion[5] = false;
//		if (regions[rgn_nghbr[idx].nb_NW].tau_veto == true && et_NW > (et_jet[idx] >> 4)) activeRegion[6] = true;
//		else activeRegion[6] = false;
//		if (regions[rgn_nghbr[idx].nb_NE].tau_veto == true && et_NE > (et_jet[idx] >> 4)) activeRegion[0] = true;
//		else activeRegion[0] = false;
//		if (regions[rgn_nghbr[idx].nb_SW].tau_veto == true && et_SW > (et_jet[idx] >> 4)) activeRegion[8] = true;
//		else activeRegion[8] = false;
//		if (regions[rgn_nghbr[idx].nb_SE].tau_veto == true && et_SE > (et_jet[idx] >> 4)) activeRegion[2] = true;
//		else activeRegion[2] = false;
//
//		rEta_jet[idx] = etapattern(activeRegion);
//		rPhi_jet[idx] = phipattern(activeRegion);
//
//		for(int iEta = 0; iEta < 3; iEta++){
//#pragma HLS UNROLL
//			bool activeStrip = false;
//			for(int iPhi = 0; iPhi < 3; iPhi++){
//#pragma HLS UNROLL
//				if(activeRegion[iEta*3+iPhi]) activeStrip = true;
//			}
//			if(activeStrip) rEta_jet[idx] |= (0x1 << iEta);
//		}
//
//		for(int iPhi = 0; iPhi < 3; iPhi++){
//#pragma HLS UNROLL
//			bool activeStrip = false;
//			for(int iEta = 0; iEta < 3; iEta++){
//#pragma HLS UNROLL
//				if(activeRegion[iEta*3+iPhi]) activeStrip = true;
//			}
//			if(activeStrip) rPhi_jet[idx] |= (0x1 << iPhi);
//		}
             
	}
}
