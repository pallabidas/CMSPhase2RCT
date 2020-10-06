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
#include "adder_tree.h"
#include "am_sort.h"
#include "PU_LUT.h"
#include "calo_out_coordinates.h"
#include "superregion.h"
#include "bitonicSort64.h"

const uint16_t NRegionsPerLink = 11; // Bits 16-31, 32-47, ..., 176-191, keeping range(15, 0) unused
const uint16_t MaxRegions = N_CH_IN * NRegionsPerLink;

 /*
  * algo_unpacked interface exposes fully unpacked input and output link data.
  * This version assumes use of 10G 8b10b links, and thus providing 192bits/BX/link.
  *
  * !!! N.B.: Do NOT use the first bytes of link_in/link_out (i.e. link_in/link_out[].range(7,0)
  * as this portion is reserved for transmission of 8b10b input/output link alignment markers.
  *
  * The remaining 184 bits are available for algorithm use.
  *
  * !!! N.B. 2: make sure to assign every bit of link_out[] data. First byte should be assigned zero.
  */

void algo_unpacked(ap_uint<192> link_in[N_CH_IN], ap_uint<192> link_out[N_CH_OUT])
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

//	for (int lnk = 0; lnk < N_CH_IN; lnk++) {
//#pragma HLS UNROLL
////  pass-through "algo"
//        link_out[lnk].range(7,0) = 0;
//        link_out[lnk].range(191,8) = link_in[lnk].range(191,8) ;
//    }
//}

	static bool first = true; //true to print 
	region_t calo_regions[NR_CALO_REG];
#pragma HLS ARRAY_PARTITION variable=calo_regions complete dim=1
	regionLoop: for(int iRegion = 0; iRegion < NR_CALO_REG; iRegion++) {
#pragma HLS UNROLL
		if(iRegion > MaxRegions) {	
			fprintf(stderr, "Too many regions - aborting");
			exit(1);
		}
		int link_idx = iRegion / NRegionsPerLink;
		int bitLo1 = ((iRegion - link_idx * NRegionsPerLink) % NRegionsPerLink + 1) * 16;
		int bitHi1 = bitLo1 + 9;
		int bitLo2 = bitHi1 + 1;
		int bitHi2 = bitLo2 + 1;
		int bitLo3 = bitHi2 + 1;
		int bitHi3 = bitLo3 + 1;
		int bitLo4 = bitHi3 + 1;
		int bitHi4 = bitLo4;
		int bitLo5 = bitHi4 + 1;
		int bitHi5 = bitLo5;
		calo_regions[iRegion].et = link_in[link_idx].range(bitHi1, bitLo1);   // 10 bits
		calo_regions[iRegion].rloc_eta = link_in[link_idx].range(bitHi2, bitLo2);   // 2 bits
		calo_regions[iRegion].rloc_phi = link_in[link_idx].range(bitHi3, bitLo3);   // 2 bits
		calo_regions[iRegion].eg_veto = link_in[link_idx].range(bitHi4, bitLo4);   // 1 bit
		calo_regions[iRegion].tau_veto = link_in[link_idx].range(bitHi5, bitLo5);   // 1 bit
		//if(first && calo_regions[iRegion].et > 0) printf("calo_regions[%d] = link_in[%d].range(%d, %d) = %d;\n", iRegion, link_idx, bitLo+9, bitLo, calo_regions[iRegion]);
		if((double)calo_regions[iRegion].et > 0) cout << "Calo region " << " ET: " << calo_regions[iRegion].et << " Eta: " << calo_regions[iRegion].rloc_eta << " Phi: " << calo_regions[iRegion].rloc_phi << " EG veto: " << calo_regions[iRegion].eg_veto << " Tau veto: " << calo_regions[iRegion].tau_veto << endl;
	}
	
////////////////////////////////////////////////////////////
	// Objets from input
	ap_uint<10> et_calo[NR_CALO_REG];
	ap_uint<10> pu_sub_et_calo[NR_CALO_REG];
	
	ap_uint<10> et_3by3_calo[NR_CALO_REG];
	
	ap_uint<10> et_jet_calo[NR_SUPER_REG];
	ap_uint<8> rIdx_jet_calo[NR_SUPER_REG];
	
	ap_uint<10> et_3by3_cntr[NR_CNTR_REG];

	ap_uint<10> et_jet_boosted[NR_SCNTR_REG];
	bitset<3>  rEta_jet_boosted[NR_SCNTR_REG];
	bitset<3>  rPhi_jet_boosted[NR_SCNTR_REG];
	ap_uint<8> rIdx_jet_boosted[NR_SCNTR_REG];
	
	ap_uint<10> nonIso_egamma_et[NR_SCNTR_REG];
	ap_uint<10> Iso_egamma_et[NR_SCNTR_REG];
	ap_uint<8> rIdx_egamma[NR_SCNTR_REG];
	
	ap_uint<10> nonIso_tau_et[NR_SCNTR_REG];
	ap_uint<10> Iso_tau_et[NR_SCNTR_REG];
	ap_uint<8> rIdx_tau[NR_SCNTR_REG];

////////////////////////////////////////////////////////////
	// Sort Objects (SO) , inputs and outputs
	ap_uint<23> so_in_jet_cr[64];
	ap_uint<23> so_out_jet_cr[64];

	ap_uint<23> so_in_jet_fwd[64];
	ap_uint<23> so_out_jet_fwd[64];

	ap_uint<23> so_in_jet_boosted[64];
	ap_uint<23> so_out_jet_boosted[64];

	ap_uint<23> so_in_tau_noniso[64];
	ap_uint<23> so_out_tau_noniso[64];

	ap_uint<23> so_in_tau_iso[64];
	ap_uint<23> so_out_tau_iso[64];

	ap_uint<23> so_in_eg_noniso[64];
	ap_uint<23> so_out_eg_noniso[64];

	ap_uint<23> so_in_eg_iso[64];
	ap_uint<23> so_out_eg_iso[64];

////////////////////////////////////////////////////////////

	ap_uint<NR_CALO_REG> tmp = 0;
	ap_uint<PUM_LEVEL_BITSIZE> pum_level;
	ap_uint<5> pum_bin;

	ap_uint<13> et_total_tmp = 0;
	ap_uint<13> et_total_ht_tmp = 0;

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

#pragma HLS ARRAY_PARTITION variable=so_in_jet_cr complete dim=0
#pragma HLS ARRAY_PARTITION variable=so_out_jet_cr complete dim=0
#pragma HLS ARRAY_PARTITION variable=so_in_jet_fwd complete dim=0
#pragma HLS ARRAY_PARTITION variable=so_out_jet_fwd complete dim=0
#pragma HLS ARRAY_PARTITION variable=so_in_jet_boosted complete dim=0
#pragma HLS ARRAY_PARTITION variable=so_out_jet_boosted complete dim=0
#pragma HLS ARRAY_PARTITION variable=so_in_tau_noniso complete dim=0
#pragma HLS ARRAY_PARTITION variable=so_out_tau_noniso complete dim=0
#pragma HLS ARRAY_PARTITION variable=so_in_tau_iso complete dim=0
#pragma HLS ARRAY_PARTITION variable=so_out_tau_iso complete dim=0
#pragma HLS ARRAY_PARTITION variable=so_in_eg_noniso complete dim=0
#pragma HLS ARRAY_PARTITION variable=so_out_eg_noniso complete dim=0
#pragma HLS ARRAY_PARTITION variable=so_in_eg_iso complete dim=0
#pragma HLS ARRAY_PARTITION variable=so_out_eg_iso complete dim=0

#pragma HLS ARRAY_PARTITION variable=et_calo complete dim=0
#pragma HLS ARRAY_PARTITION variable=pu_sub_et_calo complete dim=0
#pragma HLS ARRAY_PARTITION variable=et_3by3_calo complete dim=0
#pragma HLS ARRAY_PARTITION variable=centr_region_pu_sub complete dim=0
#pragma HLS ARRAY_PARTITION variable=et_3by3_cntr complete dim=0

#pragma HLS ARRAY_PARTITION variable=et_jet_calo complete dim=0
#pragma HLS ARRAY_PARTITION variable=rIdx_jet_calo complete dim=0

#pragma HLS ARRAY_PARTITION variable=et_jet_boosted complete dim=0
#pragma HLS ARRAY_PARTITION variable=rEta_jet_boosted complete dim=0
#pragma HLS ARRAY_PARTITION variable=rPhi_jet_boosted complete dim=0
#pragma HLS ARRAY_PARTITION variable=rIdx_jet_boosted complete dim=0

#pragma HLS ARRAY_PARTITION variable=nonIso_egamma_et complete dim=0
#pragma HLS ARRAY_PARTITION variable=Iso_egamma_et complete dim=0
#pragma HLS ARRAY_PARTITION variable=rIdx_egamma complete dim=0
#pragma HLS ARRAY_PARTITION variable=nonIso_tau_et complete dim=0
#pragma HLS ARRAY_PARTITION variable=Iso_tau_et complete dim=0
#pragma HLS ARRAY_PARTITION variable=rIdx_tau complete dim=0

////////////////////////////////////////////////////////////////////////

	//  "pum bin" calculation
	for (int i = 0; i < NR_CALO_REG; i++)
	{
#pragma HLS UNROLL
		if (calo_regions[i].et > algo_config.pum_thr)
		{
			tmp.set_bit((i), true);
		}
		else
		{
			tmp.set_bit((i), false);
		}
	}

	// count number of ones in tmp variable
	pum_level = popcount(tmp);
	pum_bin = pum_level / 26;

////////////////////////////////////////////////////////////
	// Unpack calo ET values in et_calo array
	for (int idx = 0; idx < NR_CALO_REG; idx++)
	{
#pragma HLS UNROLL
		et_calo[idx] = calo_regions[idx].et;
	}

////////////////////////////////////////////////////////////
	// Calculate pile-up subtracted ET values: pu_sub_et_calo
	pu_lut(pum_bin, et_calo, pu_sub_et_calo);

////////////////////////////////////////////////////////////
	// TODO: ener sum implementation is incomplete currently
	sums_t sums;
	sums = ener_sums(pu_sub_et_calo, 10);

////////////////////////////////////////////////////////////
	et_3by3(pu_sub_et_calo, et_3by3_calo);

	for (unsigned int phi = 0; phi < 18; phi++)
	{
#pragma HLS UNROLL

		for (int reg = 6; reg < 20; reg++)
		{
#pragma HLS UNROLL
			int idx = 26 * phi + reg;
			int cidx = 14 * phi - 6 + reg;
			et_3by3_cntr[cidx] = et_3by3_calo[idx];
			centr_region_pu_sub[cidx].et = pu_sub_et_calo[idx];
                        centr_region_pu_sub[cidx].eg_veto = calo_regions[idx].eg_veto;
                        centr_region_pu_sub[cidx].tau_veto = calo_regions[idx].tau_veto;
                        centr_region_pu_sub[cidx].rloc_eta = calo_regions[idx].rloc_eta;
                        centr_region_pu_sub[cidx].rloc_phi = calo_regions[idx].rloc_phi;
		}
	}

	// Jet algorithm
	jet(algo_config.jet_seed, pu_sub_et_calo, et_3by3_calo, et_jet_calo, rIdx_jet_calo);
	boostedjet(algo_config.jet_seed, centr_region_pu_sub, et_3by3_cntr, et_jet_boosted, rEta_jet_boosted, rPhi_jet_boosted, rIdx_jet_boosted);

	// e-gamma algorithm
	egamma(algo_config.egamma_seed, algo_config.egamma_IsoFact, centr_region_pu_sub, et_3by3_cntr, nonIso_egamma_et, Iso_egamma_et, rIdx_egamma);

	// Tau algorithm
	tau(algo_config.tau_seed, algo_config.tau_IsoFact, centr_region_pu_sub, et_3by3_cntr, nonIso_tau_et, Iso_tau_et, rIdx_tau);

////////////////////////////////////////////////////////////
	// Prepare algorithm results for sorting
	// so_in_x : sort object input
	// so_out_x: sort object output
	int bit1 = 0; int bit2 = 9; int bit3 = 10; int bit4 = 18; int bit5 = 19; int bit6 = 20; int bit7 = 21; int bit8 = 22;
	for (unsigned int phi = 0; phi < 9; phi++)
	{
#pragma HLS UNROLL
		for (int reg = 0; reg < 13; reg++)
		{
#pragma HLS UNROLL
			int idx_in = 13 * phi + reg;
			int idx_jet = rIdx_jet_calo[idx_in];

			if (reg < 3)
			{
				int idx_out = 6 * phi + reg;

				so_in_jet_fwd[idx_out].range(bit2,bit1) = et_jet_calo[idx_in];
				so_in_jet_fwd[idx_out].range(bit4, bit3) = idx_jet;
				so_in_jet_fwd[idx_out].range(bit6, bit5) = calo_regions[idx_jet].rloc_phi;
				so_in_jet_fwd[idx_out].range(bit8, bit7) = calo_regions[idx_jet].rloc_eta;
			}
			else if (reg > 9)
			{
				int idx_out = 6 * phi - 7 + reg;

				so_in_jet_fwd[idx_out].range(bit2,bit1) = et_jet_calo[idx_in];
				so_in_jet_fwd[idx_out].range(bit4, bit3) = idx_jet;
				so_in_jet_fwd[idx_out].range(bit6, bit5) = calo_regions[idx_jet].rloc_phi;
				so_in_jet_fwd[idx_out].range(bit8, bit7) = calo_regions[idx_jet].rloc_eta;
			}
			else
			{
				int idx_out = 7 * phi - 3 + reg;

				so_in_jet_cr[idx_out].range(bit2,bit1) = et_jet_calo[idx_in];
				so_in_jet_cr[idx_out].range(bit4, bit3) = idx_jet;
				so_in_jet_cr[idx_out].range(bit6, bit5) = calo_regions[idx_jet].rloc_phi;
				so_in_jet_cr[idx_out].range(bit8, bit7) = calo_regions[idx_jet].rloc_eta;
			}
		}
	}

	for (int idx = 0; idx < 63; idx++)
	{
#pragma HLS UNROLL
		int idx_egamma = rIdx_egamma[idx];
		so_in_eg_noniso[idx].range(bit2,bit1) = nonIso_egamma_et[idx];
		so_in_eg_noniso[idx].range(bit4, bit3) = idx_egamma;
		so_in_eg_noniso[idx].range(bit6, bit5) = centr_region_pu_sub[idx_egamma].rloc_phi;
		so_in_eg_noniso[idx].range(bit8, bit7) = centr_region_pu_sub[idx_egamma].rloc_eta;

		so_in_eg_iso[idx].range(bit2,bit1) = Iso_egamma_et[idx];
		so_in_eg_iso[idx].range(bit4, bit3) = idx_egamma;
		so_in_eg_iso[idx].range(bit6, bit5) = centr_region_pu_sub[idx_egamma].rloc_phi;
		so_in_eg_iso[idx].range(bit8, bit7) = centr_region_pu_sub[idx_egamma].rloc_eta;

		int idx_tau = rIdx_tau[idx];
		so_in_tau_noniso[idx].range(bit2,bit1) = nonIso_tau_et[idx];
		so_in_tau_noniso[idx].range(bit4, bit3) = idx_tau;
		so_in_tau_noniso[idx].range(bit6, bit5) = centr_region_pu_sub[idx_tau].rloc_phi;
		so_in_tau_noniso[idx].range(bit8, bit7) = centr_region_pu_sub[idx_tau].rloc_eta;

		so_in_tau_iso[idx].range(bit2,bit1) = Iso_tau_et[idx];
		so_in_tau_iso[idx].range(bit4, bit3) = idx_tau;
		so_in_tau_iso[idx].range(bit6, bit5) = centr_region_pu_sub[idx_tau].rloc_phi;
		so_in_tau_iso[idx].range(bit8, bit7) = centr_region_pu_sub[idx_tau].rloc_eta;

		int idx_boosted = rIdx_jet_boosted[idx];
		so_in_jet_boosted[idx].range(bit2,bit1) = et_jet_boosted[idx];
		so_in_jet_boosted[idx].range(bit4, bit3) = idx_boosted;
		so_in_jet_boosted[idx].range(bit6, bit5) = centr_region_pu_sub[idx_boosted].rloc_phi;
		so_in_jet_boosted[idx].range(bit8, bit7) = centr_region_pu_sub[idx_boosted].rloc_eta;
	}

	for (int idx = 54; idx < 64; idx++)
	{
#pragma HLS UNROLL

		//so_in_jet_fwd[idx].range(bit2, bit1) = 0;
		//so_in_jet_fwd[idx].range(bit4, bit3) = 0;
		//so_in_jet_fwd[idx].range(bit6, bit5) = 0;
		//so_in_jet_fwd[idx].range(bit8, bit7) = 0;
		so_in_jet_fwd[idx] = 0;
	}


	for (int idx = 63; idx < 64; idx++)
	{
#pragma HLS UNROLL
		so_in_jet_cr[idx] = 0;
		so_in_eg_noniso[idx] = 0;
		so_in_eg_iso[idx] = 0;
		so_in_tau_noniso[idx] = 0;
		so_in_tau_iso[idx] = 0;
		so_in_jet_boosted[idx] = 0;
		//so_in_jet_cr[idx].range(bit2, bit1) = 0;
		//so_in_jet_cr[idx].range(bit4, bit3) = 0;
		//so_in_jet_cr[idx].range(bit6, bit5) = 0;
		//so_in_jet_cr[idx].range(bit8, bit7) = 0;
		//
		//so_in_eg_noniso[idx].range(bit2, bit1) = 0;
		//so_in_eg_noniso[idx].range(bit4, bit3) = 0;
		//so_in_eg_noniso[idx].range(bit6, bit5) = 0;
		//so_in_eg_noniso[idx].range(bit8, bit7) = 0;
		//
		//so_in_eg_iso[idx].range(bit2, bit1) = 0;
		//so_in_eg_iso[idx].range(bit4, bit3) = 0;
		//so_in_eg_iso[idx].range(bit6, bit5) = 0;
		//so_in_eg_iso[idx].range(bit8, bit7) = 0;
		//
		//so_in_tau_noniso[idx].range(bit2, bit1) = 0;
		//so_in_tau_noniso[idx].range(bit4, bit3) = 0;
		//so_in_tau_noniso[idx].range(bit6, bit5) = 0;
		//so_in_tau_noniso[idx].range(bit8, bit7) = 0;
		//
		//so_in_tau_iso[idx].range(bit2, bit1) = 0;
		//so_in_tau_iso[idx].range(bit4, bit3) = 0;
		//so_in_tau_iso[idx].range(bit6, bit5) = 0;
		//so_in_tau_iso[idx].range(bit8, bit7) = 0;

		//so_in_jet_boosted[idx].range(bit2, bit1) = 0;
		//so_in_jet_boosted[idx].range(bit4, bit3) = 0;
		//so_in_jet_boosted[idx].range(bit6, bit5) = 0;
		//so_in_jet_boosted[idx].range(bit8, bit7) = 0;
	}	

////////////////////////////////////////////////////////////
	// Do the actual sorting
	bitonicSort64(so_in_jet_cr, so_out_jet_cr);
	bitonicSort64(so_in_jet_fwd, so_out_jet_fwd);
	bitonicSort64(so_in_jet_boosted, so_out_jet_boosted);

	bitonicSort64(so_in_eg_noniso, so_out_eg_noniso);
	bitonicSort64(so_in_eg_iso, so_out_eg_iso);

	bitonicSort64(so_in_tau_noniso, so_out_tau_noniso);
	bitonicSort64(so_in_tau_iso, so_out_tau_iso);


////////////////////////////////////////////////////////////
	// Assign the top 8 candidates to algorithnm outputs
	for (int idx = 63; idx > 55; idx--)
	{
#pragma HLS UNROLL
		ap_uint<1> side;
		ap_uint<9> idx_srt;
		ap_uint<6> ieta;
		int idx_out = 63 - idx;


		{ // Central Jets
			int bLo1 = 16;
			int bHi1 = bLo1 + 10; // 11 bits
			tmp_link_out[idx_out].range(bHi1, bLo1) = so_out_jet_cr[idx].range(bit2, bit1);
			idx_srt = so_out_jet_cr[idx].range(bit4, bit3);

			side = calo_coor_full[idx_srt].side;
			int bLo2 = bHi1 + 1;
			int bHi2 = bLo2; // 1 bit
			tmp_link_out[idx_out].range(bHi2, bLo2) = side;

			int bLo3 = bHi2 + 1;
			int bHi3 = bLo3 + 6; // 7 bits
			tmp_link_out[idx_out].range(bHi3, bLo3) = calo_coor_full[idx_srt].iphi + so_out_jet_cr[idx].range(bit6, bit5);

			int bLo4 = bHi3 + 1;
			int bHi4 = bLo4 + 5; // 6 bits
			if (side == 1)
				tmp_link_out[idx_out].range(bHi4, bLo4) = calo_coor_full[idx_srt].ieta - so_out_jet_cr[idx].range(bit8, bit7);
			else
				tmp_link_out[idx_out].range(bHi4, bLo4) = calo_coor_full[idx_srt].ieta + so_out_jet_cr[idx].range(bit8, bit7);

			if((double)tmp_link_out[idx_out].range(bHi1, bLo1) > 0) cout << "Jet CR " << idx_out << " ET: " << dec << tmp_link_out[idx_out].range(bHi1, bLo1) << " Side: " << tmp_link_out[idx_out].range(bHi2, bLo2) << " iPhi: " << tmp_link_out[idx_out].range(bHi3, bLo3)  << " iEta: " << tmp_link_out[idx_out].range(bHi4, bLo4) << endl ;
		}

		{ // Forward Jets
			int bLo5 = 16 + 25 + 1;
			int bHi5 = bLo5 + 10;
			tmp_link_out[idx_out].range(bHi5, bLo5) = so_out_jet_fwd[idx].range(bit2, bit1);
			idx_srt = so_out_jet_fwd[idx].range(bit4, bit3);

			side = calo_coor_full[idx_srt].side;
			int bLo6 = bHi5 + 1;
			int bHi6 = bLo6;
			tmp_link_out[idx_out].range(bHi6, bLo6) = side;

			int bLo7 = bHi6 + 1;
			int bHi7 = bLo7 + 6; 
			tmp_link_out[idx_out].range(bHi7, bLo7) = calo_coor_full[idx_srt].iphi + so_out_jet_fwd[idx].range(bit6, bit5);

			int bLo8 = bHi7 + 1;
			int bHi8 = bLo8 + 5;
			if (side == 1)
				tmp_link_out[idx_out].range(bHi8, bLo8) = calo_coor_full[idx_srt].ieta - so_out_jet_fwd[idx].range(bit8, bit7);
			else
				tmp_link_out[idx_out].range(bHi8, bLo8) = calo_coor_full[idx_srt].ieta + so_out_jet_fwd[idx].range(bit8, bit7);

			if((double)tmp_link_out[idx_out].range(bHi5, bLo5) > 0) cout << "Jet FWD " << idx_out << " ET: " << dec << tmp_link_out[idx_out].range(bHi5, bLo5) << " Side: " << tmp_link_out[idx_out].range(bHi6, bLo6) << " iPhi: " << tmp_link_out[idx_out].range(bHi7, bLo7)  << " iEta: " << tmp_link_out[idx_out].range(bHi8, bLo8) << endl ;
		}
		
		{ // Boosted jets
			int bLo9 = 16 + 25 + 25 + 1;
			int bHi9 = bLo9 + 10;
			tmp_link_out[idx_out].range(bHi9, bLo9) = so_out_jet_boosted[idx].range(bit2, bit1);
			idx_srt = so_out_jet_boosted[idx].range(bit4, bit3);
			
			side = calo_coor[idx_srt].side;
			int bLo10 = bHi9 + 1;
			int bHi10 = bLo10;
			tmp_link_out[idx_out].range(bHi10, bLo10) = side;
	
			int bLo11 = bHi10 + 1;
			int bHi11 = bLo11 + 6;
			tmp_link_out[idx_out].range(bHi11, bLo11) = calo_coor[idx_srt].iphi + so_out_jet_boosted[idx].range(bit6, bit5);

			int bLo12 = bHi11 + 1;
			int bHi12 = bLo12 + 5;
			if (side == 1)
				tmp_link_out[idx_out].range(bHi12, bLo12) = calo_coor[idx_srt].ieta - so_out_jet_boosted[idx].range(bit8, bit7);
			else
				tmp_link_out[idx_out].range(bHi12, bLo12) = calo_coor[idx_srt].ieta + so_out_jet_boosted[idx].range(bit8, bit7);

			int sidx = central_super_region_idx[idx_srt];
			int bLo13 = bHi12 + 1;
			int bHi13 = bLo13 + 2;
			tmp_link_out[idx_out].range(bHi13, bLo13) = rPhi_jet_boosted[sidx].to_ulong();

			int bLo14 = bHi13 + 1;
			int bHi14 = bLo14 + 2;
			tmp_link_out[idx_out].range(bHi14, bLo14) = rEta_jet_boosted[sidx].to_ulong();

			if((double)tmp_link_out[idx_out].range(bHi9, bLo9) > 0) cout << "Jet Boosted " << idx_out << " ET: " << dec << tmp_link_out[idx_out].range(bHi9, bLo9) << " Side: " << tmp_link_out[idx_out].range(bHi10, bLo10) << " iPhi: " << tmp_link_out[idx_out].range(bHi11, bLo11)  << " iEta: " << tmp_link_out[idx_out].range(bHi12, bLo12) << " rPhi: " << tmp_link_out[idx_out].range(bHi13, bLo13).to_ulong() << " rEta: " << tmp_link_out[idx_out].range(bHi14, bLo14).to_ulong() << endl;
		}

		{ // EG NonIso
			int bLo15 = 16 + 25 + 25 + 31 + 1;
			int bHi15 = bLo15 + 8;
			tmp_link_out[idx_out].range(bHi15, bLo15) = so_out_eg_noniso[idx].range(bit2, bit1);
			idx_srt = so_out_eg_noniso[idx].range(bit4, bit3);

			side = calo_coor[idx_srt].side;
			int bLo16 = bHi15 + 1;
			int bHi16 = bLo16;
			tmp_link_out[idx_out].range(bHi16, bLo16) = side;

			int bLo17 = bHi16 + 1;
			int bHi17 = bLo17 + 6;
			tmp_link_out[idx_out].range(bHi17, bLo17) = calo_coor[idx_srt].iphi + so_out_eg_noniso[idx].range(bit6, bit5);

			int bLo18 = bHi17 + 1;
			int bHi18 = bLo18 + 5;
			if (side == 1)
				tmp_link_out[idx_out].range(bHi18, bLo18) = calo_coor[idx_srt].ieta - so_out_eg_noniso[idx].range(bit8, bit7);
			else
				tmp_link_out[idx_out].range(bHi18, bLo18) = calo_coor[idx_srt].ieta + so_out_eg_noniso[idx].range(bit8, bit7);

			if((double)tmp_link_out[idx_out].range(bHi15, bLo15) > 0) cout << "EG NonIso " << idx_out << " ET: " << dec << tmp_link_out[idx_out].range(bHi15, bLo15) << " Side: " << tmp_link_out[idx_out].range(bHi16, bLo16) << " iPhi: " << tmp_link_out[idx_out].range(bHi17, bLo17)  << " iEta: " << tmp_link_out[idx_out].range(bHi18, bLo18) << endl ;
		}

		{ //EG ISO
			int bLo19 = 16 + 25 + 25 + 31 + 23 + 1;
			int bHi19 = bLo19 + 8;
			tmp_link_out[idx_out].range(bHi19, bLo19) = so_out_eg_iso[idx].range(bit2, bit1);
			idx_srt = so_out_eg_iso[idx].range(bit4, bit3);

			side = calo_coor[idx_srt].side;
			int bLo20 = bHi19 + 1;
			int bHi20 = bLo20;
			tmp_link_out[idx_out].range(bHi20, bLo20) = side;

			int bLo21 = bHi20 + 1;
			int bHi21 = bLo21 + 6;
			tmp_link_out[idx_out].range(bHi21, bLo21) = calo_coor[idx_srt].iphi + so_out_eg_iso[idx].range(bit6, bit5);

			int bLo22 = bHi21 + 1;
			int bHi22 = bLo22 + 5;
			if (side == 1)
				tmp_link_out[idx_out].range(bHi22, bLo22) = calo_coor[idx_srt].ieta - so_out_eg_iso[idx].range(bit8, bit7);
			else
				tmp_link_out[idx_out].range(bHi22, bLo22) = calo_coor[idx_srt].ieta + so_out_eg_iso[idx].range(bit8, bit7);

			if((double)tmp_link_out[idx_out].range(bHi19, bLo19) > 0) cout << "EG Iso " << idx_out << " ET: " << dec << tmp_link_out[idx_out].range(bHi19, bLo19) << " Side: " << tmp_link_out[idx_out].range(bHi20, bLo20) << " iPhi: " << tmp_link_out[idx_out].range(bHi21, bLo21)  << " iEta: " << tmp_link_out[idx_out].range(bHi22, bLo22) << endl ;
		}

		{ // Tau NonIso
			int bLo23 = 16 + 25 + 25 + 31 + 23 + 23 + 1;
			int bHi23 = bLo23 + 8;
			tmp_link_out[idx_out].range(bHi23, bLo23) = so_out_tau_noniso[idx].range(bit2, bit1);
			idx_srt = so_out_tau_noniso[idx].range(bit4, bit3);

			side = calo_coor[idx_srt].side;
			int bLo24 = bHi23 + 1;
			int bHi24 = bLo24;
			tmp_link_out[idx_out].range(bHi24, bLo24) = side;

			int bLo25 = bHi24 + 1;
			int bHi25 = bLo25 + 6;
			tmp_link_out[idx_out].range(bHi25, bLo25) = calo_coor[idx_srt].iphi + so_out_tau_noniso[idx].range(bit6, bit5);

			int bLo26 = bHi25 + 1;
			int bHi26 = bLo26 + 5;
			if (side == 1)
				tmp_link_out[idx_out].range(bHi26, bLo26) = calo_coor[idx_srt].ieta - so_out_tau_noniso[idx].range(bit8, bit7);
			else
				tmp_link_out[idx_out].range(bHi26, bLo26) = calo_coor[idx_srt].ieta + so_out_tau_noniso[idx].range(bit8, bit7);

			if((double)tmp_link_out[idx_out].range(bHi23, bLo23) > 0) cout << "Tau NonIso " << idx_out << " ET: " << dec << tmp_link_out[idx_out].range(bHi23, bLo23) << " Side: " << tmp_link_out[idx_out].range(bHi24, bLo24) << " iPhi: " << tmp_link_out[idx_out].range(bHi25, bLo25)  << " iEta: " << tmp_link_out[idx_out].range(bHi26, bLo26) << endl ;
		}

		{ //Tau ISO
			int bLo27 = 16 + 25 + 25 + 31 + 23 + 23 + 23 + 1;
			int bHi27 = bLo27 + 8;
			tmp_link_out[idx_out].range(bHi27, bLo27) = so_out_tau_iso[idx].range(bit2, bit1);
			idx_srt = so_out_tau_iso[idx].range(bit4, bit3);

			side = calo_coor[idx_srt].side;
			int bLo28 = bHi27 + 1;
			int bHi28 = bLo28;
			tmp_link_out[idx_out].range(bHi28, bLo28) = side;

			int bLo29 = bHi28 + 1;
			int bHi29 = bLo29 + 6;
			tmp_link_out[idx_out].range(bHi29, bLo29) = calo_coor[idx_srt].iphi + so_out_tau_iso[idx].range(bit6, bit5);

			int bLo30 = bHi29 + 1;
			int bHi30 = bLo30 + 5;
			if (side == 1)
				tmp_link_out[idx_out].range(bHi30, bLo30) = calo_coor[idx_srt].ieta - so_out_tau_iso[idx].range(bit8, bit7);
			else
				tmp_link_out[idx_out].range(bHi30, bLo30) = calo_coor[idx_srt].ieta + so_out_tau_iso[idx].range(bit8, bit7);

			if((double)tmp_link_out[idx_out].range(bHi27, bLo27) > 0) cout << "Tau Iso " << idx_out << " ET: " << dec << tmp_link_out[idx_out].range(bHi27, bLo27) << " Side: " << tmp_link_out[idx_out].range(bHi28, bLo28) << " iPhi: " << tmp_link_out[idx_out].range(bHi29, bLo29)  << " iEta: " << tmp_link_out[idx_out].range(bHi30, bLo30) << endl ;
		}
	}

	for(int i = 0; i < N_CH_OUT; i++){
#pragma HLS unroll
		link_out[i] = tmp_link_out[i];
	}
}

////////////////////////////////////////////////////////////
// count number of ones in bitString
ap_uint<8> popcount(ap_uint<NR_CALO_REG> bitString)
{
        ap_uint<9> popcnt = 0;

        for (ap_uint<9> b = 0; b < NR_CALO_REG; b++)
        {
#pragma HLS unroll
                popcnt += ((bitString >> b) & 1);
        }
        return popcnt;
}

