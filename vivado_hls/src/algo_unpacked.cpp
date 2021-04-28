#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <unistd.h>

#include <iostream>
#include <fstream>
#include <iomanip>
#include <string>

using namespace std;

#include "algo_unpacked.h"
#include "UCTSummaryCard.hpp"
#include "PU_LUT.h"
#include "calo_out_coordinates.h"
#include "superregion.h"
#include "bitonicSort64.h"

const uint16_t NRegionsPerLink = 7; // Bits 16-31, 32-47, ..., 96-112, keeping range(15, 0) unused
const uint16_t MaxRegions = N_CH_IN * NRegionsPerLink;

 /*
  * algo_unpacked interface exposes fully unpacked input and output link data.
  * This version assumes use of 8G 8b10b links, and thus providing 192bits/BX/link.
  *
  * !!! N.B.: Do NOT use the first bytes of link_in/link_out (i.e. link_in/link_out[].range(7,0)
  * as this portion is reserved for transmission of 8b10b input/output link alignment markers.
  *
  * The remaining 184 bits are available for algorithm use.
  *
  * !!! N.B. 2: make sure to assign every bit of link_out[] data. First byte should be assigned zero.
  */

void algo_unpacked(ap_uint<112> link_in[N_CH_IN], ap_uint<192> link_out[N_CH_OUT])
{

// !!! Retain these 4 #pragma directives below in your algo_unpacked implementation !!!
#pragma HLS ARRAY_PARTITION variable=link_in complete dim=0
#pragma HLS ARRAY_PARTITION variable=link_out complete dim=0
#pragma HLS PIPELINE II=3
#pragma HLS INTERFACE ap_ctrl_hs port=return

	ap_uint<192> tmp_link_out[N_CH_OUT];
#pragma HLS ARRAY_PARTITION variable=tmp_link_out    complete dim=0
	for (int idx = 0; idx < N_CH_OUT; idx++){
#pragma HLS UNROLL
		tmp_link_out[idx]         = 0;
	}   

	// null algo specific pragma: avoid fully combinatorial algo by specifying min latency
	// otherwise algorithm clock input (ap_clk) gets optimized away
#pragma HLS latency min=3

	static bool first = true; //true to print 
	region_t centr_region[NR_CNTR_REG];
#pragma HLS ARRAY_PARTITION variable=centr_region complete dim=1
	regionLoop: for(int iRegion = 0; iRegion < NR_CNTR_REG; iRegion++) {
#pragma HLS UNROLL
		if(iRegion > MaxRegions) {	
			fprintf(stderr, "Too many regions - aborting");
			exit(1);
		}
		int link_idx = iRegion / NRegionsPerLink;
		int bitLo1 = ((iRegion - link_idx * NRegionsPerLink) % NRegionsPerLink) * 16;
		int bitHi1 = bitLo1 + 9;
		int bitLo2 = bitHi1 + 1;
		int bitHi2 = bitLo2 + 1;
		int bitLo3 = bitHi2 + 1;
		int bitHi3 = bitLo3 + 1;
		int bitLo4 = bitHi3 + 1;
		int bitHi4 = bitLo4;
		int bitLo5 = bitHi4 + 1;
		int bitHi5 = bitLo5;
		centr_region[iRegion].et = link_in[link_idx].range(bitHi1, bitLo1);   // 10 bits
		centr_region[iRegion].rloc_eta = link_in[link_idx].range(bitHi2, bitLo2);   // 2 bits
		centr_region[iRegion].rloc_phi = link_in[link_idx].range(bitHi3, bitLo3);   // 2 bits
		centr_region[iRegion].eg_veto = link_in[link_idx].range(bitHi4, bitLo4);   // 1 bit
		centr_region[iRegion].tau_veto = link_in[link_idx].range(bitHi5, bitLo5);   // 1 bit
		if((double)centr_region[iRegion].et > 0) cout << "Calo region " << " ET: " << centr_region[iRegion].et << " Eta: " << centr_region[iRegion].rloc_eta << " Phi: " << centr_region[iRegion].rloc_phi << " EG veto: " << centr_region[iRegion].eg_veto << " Tau veto: " << centr_region[iRegion].tau_veto << endl;
	}
      
////////////////////////////////////////////////////////////
	// Objets from input
	ap_uint<10> et_calo[NR_CNTR_REG];
	ap_uint<10> pu_sub_et_calo[NR_CNTR_REG];	
	ap_uint<10> et_3by3_calo[NR_CNTR_REG];
	ap_uint<10> et_3by3_cntr[NR_CNTR_REG];

	ap_uint<10> et_jet_boosted[NR_SCNTR_REG];
	bitset<3>  rEta_jet_boosted[NR_SCNTR_REG];
	bitset<3>  rPhi_jet_boosted[NR_SCNTR_REG];
	ap_uint<9> rIdx_boostedjet[NR_SCNTR_REG];
	
	ap_uint<32> so_in_jet_boosted[64];
	ap_uint<32> so_out_jet_boosted[64];

	ap_uint<NR_CNTR_REG> tmp = 0;
	ap_uint<PUM_LEVEL_BITSIZE> pum_level;
	ap_uint<5> pum_bin;

	region_t centr_region_pu_sub[NR_CNTR_REG];

///////////////////////////////////////////////////////////

	algo_config_t algo_config;
	
	algo_config.egamma_IsoFact =  0.3;
	algo_config.egamma_seed = 5;
	algo_config.jet_seed = 10;
	algo_config.pum_thr = 0;
	algo_config.tau_IsoFact =  0.3;
	algo_config.tau_seed = 10;

///////////////////////////////////////////////////////////

#pragma HLS INTERFACE ap_none port=algo_config

#pragma HLS ARRAY_RESHAPE variable=centr_region_pu_sub complete dim=1

#pragma HLS ARRAY_RESHAPE variable=so_in_jet_boosted complete dim=0
#pragma HLS ARRAY_RESHAPE variable=so_out_jet_boosted complete dim=0

#pragma HLS ARRAY_RESHAPE variable=et_calo complete dim=0
#pragma HLS ARRAY_RESHAPE variable=pu_sub_et_calo complete dim=0

#pragma HLS ARRAY_RESHAPE variable=et_jet_boosted complete dim=0
#pragma HLS ARRAY_RESHAPE variable=rEta_jet_boosted complete dim=0
#pragma HLS ARRAY_RESHAPE variable=rPhi_jet_boosted complete dim=0
#pragma HLS ARRAY_RESHAPE variable=rIdx_boostedjet complete dim=0

////////////////////////////////////////////////////////////////////////
	//  "pum bin" calculation
	for (int i = 0; i < NR_CNTR_REG; i++)
	{
#pragma HLS UNROLL
		if (centr_region[i].et > algo_config.pum_thr)
		{
			tmp.set_bit((i), true);
		}
		else
		{
			tmp.set_bit((i), false);
		}
	}

	// Count number of ones in tmp variable
	pum_level = popcount(tmp);
	pum_bin = pum_level / 14;

////////////////////////////////////////////////////////////
	// Unpack calo ET values in et_calo array
	for (int idx = 0; idx < NR_CNTR_REG; idx++)
	{
#pragma HLS UNROLL
		et_calo[idx] = centr_region[idx].et;
	}

////////////////////////////////////////////////////////////
	// Calculate pile-up subtracted ET values: pu_sub_et_calo
	pu_lut_cntr(pum_bin, et_calo, pu_sub_et_calo);

////////////////////////////////////////////////////////////
	et_3by3(pu_sub_et_calo, et_3by3_calo);

	for (int idx = 0; idx < NR_CNTR_REG; idx++)
	{
#pragma HLS UNROLL
		centr_region_pu_sub[idx].et = pu_sub_et_calo[idx];
		centr_region_pu_sub[idx].eg_veto = centr_region[idx].eg_veto;
		centr_region_pu_sub[idx].tau_veto = centr_region[idx].tau_veto;
		centr_region_pu_sub[idx].rloc_eta = centr_region[idx].rloc_eta;
		centr_region_pu_sub[idx].rloc_phi = centr_region[idx].rloc_phi;
	}

////////////////////////////////////////////////////////////
	// Prepare algorithm results
	boostedjet(algo_config.jet_seed, centr_region_pu_sub, et_3by3_calo, et_jet_boosted, rEta_jet_boosted, rPhi_jet_boosted, rIdx_boostedjet);

	for (int idx = 0; idx < NR_SCNTR_REG; idx++)
	{
#pragma HLS UNROLL
		ap_uint<9> idx_jet_in = rIdx_boostedjet[idx];
		so_in_jet_boosted[idx].range(9, 0) = et_jet_boosted[idx];
		so_in_jet_boosted[idx].range(18, 10) = idx_jet_in;
		so_in_jet_boosted[idx].range(25, 19) = centr_region[idx_jet_in].rloc_phi;
		so_in_jet_boosted[idx].range(31, 26) = centr_region[idx_jet_in].rloc_eta;
	}
	so_in_jet_boosted[63] = 0;

	// Sorting objects
	bitonicSort64(so_in_jet_boosted, so_out_jet_boosted);

	// Assign the algorithnm outputs
	for (int idx = 0; idx < 12; idx++)
	{
#pragma HLS UNROLL
		ap_uint<1> side;
		ap_uint<9> idx_srt;
		
		{ // Boosted jets
			int bLo9 = 0;
			int bHi9 = bLo9 + 9;
			tmp_link_out[idx].range(bHi9, bLo9) = so_out_jet_boosted[idx].range(9, 0);
			idx_srt = so_out_jet_boosted[idx].range(18, 10);
			
			side = calo_coor[idx_srt].side;
			int bLo10 = bHi9 + 1;
			int bHi10 = bLo10;
			tmp_link_out[idx].range(bHi10, bLo10) = side;
	
			int bLo11 = bHi10 + 1;
			int bHi11 = bLo11 + 6;
			tmp_link_out[idx].range(bHi11, bLo11) = calo_coor[idx_srt].iphi + so_out_jet_boosted[idx].range(25, 19);

			int bLo12 = bHi11 + 1;
			int bHi12 = bLo12 + 5;
			if (side == 1)
				tmp_link_out[idx].range(bHi12, bLo12) = calo_coor[idx_srt].ieta - so_out_jet_boosted[idx].range(31, 26);
			else
				tmp_link_out[idx].range(bHi12, bLo12) = calo_coor[idx_srt].ieta + so_out_jet_boosted[idx].range(31, 26);

			if((double)tmp_link_out[idx].range(bHi9, bLo9) > 0) cout << "Boosted jet :" << idx << " ET: " << dec << tmp_link_out[idx].range(bHi9, bLo9) << " Side: " << tmp_link_out[idx].range(bHi10, bLo10) << " iPhi: " << tmp_link_out[idx].range(bHi11, bLo11)  << " iEta: " << tmp_link_out[idx].range(bHi12, bLo12) << endl;
		}
	}

	for(int i = 0; i < N_CH_OUT; i++){
#pragma HLS unroll
		link_out[i] = tmp_link_out[i];
	}
}

////////////////////////////////////////////////////////////
// count number of ones in bitString
ap_uint<8> popcount(ap_uint<NR_CNTR_REG> bitString)
{
        ap_uint<9> popcnt = 0;

        for (ap_uint<9> b = 0; b < NR_CNTR_REG; b++)
        {
#pragma HLS unroll
                popcnt += ((bitString >> b) & 1);
        }
        return popcnt;
}

