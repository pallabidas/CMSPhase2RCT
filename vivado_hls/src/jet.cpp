#include "ap_int.h"

#include "region_neighbors.h"

void jet(ap_uint<10> jet_seed,             // input
			  ap_uint<10> et_rgn [NR_CALO_REG], // input 26x18
			  ap_uint<10> et_3by3[NR_CALO_REG], // input 26x18
			  ap_uint<10> et_jet [NR_CALO_REG]) // *output* 26x18
{

	bool jet_veto[NR_CALO_REG];

#pragma HLS INTERFACE ap_none port=jet_seed

#pragma HLS PIPELINE II=6 // target clk freq = 250 MHz

#pragma HLS ARRAY_RESHAPE  variable=et_rgn    complete  dim=1
#pragma HLS ARRAY_RESHAPE  variable=et_3by3   complete  dim=1
#pragma HLS ARRAY_RESHAPE  variable=et_jet    complete  dim=1
#pragma HLS ARRAY_RESHAPE  variable=jet_veto  complete  dim=1

	loop_rgn_et: for (int idx = 0; idx < NR_CALO_REG; idx++)
	{
#pragma HLS UNROLL
		ap_uint<10> et_C, et_N, et_E, et_W, et_S,
					et_NE, et_NW, et_SE, et_SW;

		et_C = et_rgn[rgn_nghbr[idx].nb_C]; // center
		et_N = et_rgn[rgn_nghbr[idx].nb_N]; // north
		et_S = et_rgn[rgn_nghbr[idx].nb_S]; // south

		// do we have a border on the east?
		if (rgn_nghbr[idx].nb_E != -1) {
			et_E  = et_rgn[rgn_nghbr[idx].nb_E];  // east
			et_NE = et_rgn[rgn_nghbr[idx].nb_NE]; // north east
			et_SE = et_rgn[rgn_nghbr[idx].nb_SE]; // south east
		}
		else {
			et_E  = 0;
			et_NE = 0;
			et_SE = 0;
		}

		// do we have a border on the west?
		if (rgn_nghbr[idx].nb_W != -1) {
			et_W  = et_rgn[rgn_nghbr[idx].nb_W];  // west
			et_NW = et_rgn[rgn_nghbr[idx].nb_NW]; // north west
			et_SW = et_rgn[rgn_nghbr[idx].nb_SW]; // south west
		}
		else {
			et_W  = 0;
			et_NW = 0;
			et_SW = 0;
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

        // assign et_jet
		if (jet_veto[idx] == false)
			et_jet[idx] = et_3by3[idx];
		else
			et_jet[idx] = 0;
	}
}