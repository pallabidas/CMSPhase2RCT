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
	
	ap_uint<10> et_jet_calo[NR_CALO_REG];
	
	ap_uint<10> et_3by3_cntr[252];
	ap_uint<10> et_jet_cntr[216];
	ap_uint<10> et_jet_fwd[252];
	ap_uint<10> et_jet_boosted[252];
	bitset<3>  rEta_jet_boosted[252];
	bitset<3>  rPhi_jet_boosted[252];
	
	ap_uint<10> nonIso_egamma_et[252];
	ap_uint<10> Iso_egamma_et[252];
	
	ap_uint<10> nonIso_tau_et[252];
	ap_uint<10> Iso_tau_et[252];

////////////////////////////////////////////////////////////
	// Sort Objects (SO) , inputs and outputs
	t_so so_in_jet_cr[256];
	t_so so_out_jet_cr[8];

	t_so so_in_jet_fwd[256];
	t_so so_out_jet_fwd[8];

	t_so so_in_jet_boosted[256];
	t_so so_out_jet_boosted[8];

	t_so so_in_tau_noniso[256];
	t_so so_out_tau_noniso[8];

	t_so so_in_tau_iso[256];
	t_so so_out_tau_iso[8];

	t_so so_in_eg_noniso[256];
	t_so so_out_eg_noniso[8];

	t_so so_in_eg_iso[256];
	t_so so_out_eg_iso[8];

////////////////////////////////////////////////////////////

	ap_uint<NR_CALO_REG> tmp = 0;
	ap_uint<PUM_LEVEL_BITSIZE> pum_level;
	ap_uint<5> pum_bin;

	ap_uint<13> et_total_tmp = 0;
	ap_uint<13> et_total_ht_tmp = 0;

	region_t centr_region[NR_CNTR_REG];
	region_t fwd_region[NR_FWD_REG];
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

#pragma HLS ARRAY_RESHAPE variable=centr_region complete dim=1
#pragma HLS ARRAY_RESHAPE variable=fwd_region complete dim=1
#pragma HLS ARRAY_RESHAPE variable=centr_region_pu_sub complete dim=1

#pragma HLS ARRAY_RESHAPE variable=so_in_jet_cr complete dim=0
#pragma HLS ARRAY_RESHAPE variable=so_out_jet_cr complete dim=0
#pragma HLS ARRAY_RESHAPE variable=so_in_jet_fwd complete dim=0
#pragma HLS ARRAY_RESHAPE variable=so_out_jet_fwd complete dim=0
#pragma HLS ARRAY_RESHAPE variable=so_in_jet_boosted complete dim=0
#pragma HLS ARRAY_RESHAPE variable=so_out_jet_boosted complete dim=0
#pragma HLS ARRAY_RESHAPE variable=so_in_tau_noniso complete dim=0
#pragma HLS ARRAY_RESHAPE variable=so_out_tau_noniso complete dim=0
#pragma HLS ARRAY_RESHAPE variable=so_in_tau_iso complete dim=0
#pragma HLS ARRAY_RESHAPE variable=so_out_tau_iso complete dim=0
#pragma HLS ARRAY_RESHAPE variable=so_in_eg_noniso complete dim=0
#pragma HLS ARRAY_RESHAPE variable=so_out_eg_noniso complete dim=0
#pragma HLS ARRAY_RESHAPE variable=so_in_eg_iso complete dim=0
#pragma HLS ARRAY_RESHAPE variable=so_out_eg_iso complete dim=0

#pragma HLS ARRAY_RESHAPE variable=et_calo complete dim=0
#pragma HLS ARRAY_RESHAPE variable=pu_sub_et_calo complete dim=0

#pragma HLS ARRAY_RESHAPE variable=et_jet_cntr complete dim=0
#pragma HLS ARRAY_RESHAPE variable=et_jet_fwd complete dim=0

#pragma HLS ARRAY_RESHAPE variable=et_jet_calo complete dim=0

#pragma HLS ARRAY_RESHAPE variable=et_jet_boosted complete dim=0
#pragma HLS ARRAY_RESHAPE variable=rEta_jet_boosted complete dim=0
#pragma HLS ARRAY_RESHAPE variable=rPhi_jet_boosted complete dim=0

#pragma HLS ARRAY_RESHAPE variable=nonIso_egamma_et complete dim=0
#pragma HLS ARRAY_RESHAPE variable=Iso_egamma_et complete dim=0
#pragma HLS ARRAY_RESHAPE variable=nonIso_tau_et complete dim=0
#pragma HLS ARRAY_RESHAPE variable=Iso_tau_et complete dim=0

////////////////////////////////////////////////////////////////////////

	// unpack central and forward regions into combined "calo regions"
	for (unsigned int phi = 0; phi < 18; phi++)
	{
#pragma HLS UNROLL
		// 6 FWD-, 7 CNTR-, 7 CNTR+ and 6 FWD+ calo regions, 26 in total
		for (int reg = 0; reg < 26; reg++)
		{
#pragma HLS UNROLL
			if (reg <= 5) // FWD-
			{
				fwd_region[12 * phi + reg] = calo_regions[26 * phi + reg];
			}
			else if (reg >= 20) // FWD+
			{
				fwd_region[12 * phi - 14 + reg] = calo_regions[26 * phi + reg];
			}
			else // CNTR- and CNTR+
			{
				centr_region[14 * phi - 6 + reg] = calo_regions[26 * phi + reg];
			}
		}
	}

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

	for (int idx = 0; idx < 252; idx++)
	{
#pragma HLS UNROLL

		centr_region_pu_sub[idx].eg_veto = centr_region[idx].eg_veto;
		centr_region_pu_sub[idx].tau_veto = centr_region[idx].tau_veto;
		centr_region_pu_sub[idx].rloc_eta = centr_region[idx].rloc_eta;
		centr_region_pu_sub[idx].rloc_phi = centr_region[idx].rloc_phi;
	}

	for (unsigned int phi = 0; phi < 18; phi++)
	{
#pragma HLS UNROLL

		for (int reg = 0; reg < 26; reg++)
		{
#pragma HLS UNROLL

			if (reg > 5 && reg < 20)
			{
				et_3by3_cntr[phi * 14 + reg - 6] = et_3by3_calo[phi * 26 + reg];
				centr_region_pu_sub[phi * 14 + reg - 6].et = pu_sub_et_calo[phi	* 26 + reg];
			}
		}
	}


	// Jet algorithm
	jet(algo_config.jet_seed, pu_sub_et_calo, et_3by3_calo, et_jet_calo);
	boostedjet(algo_config.jet_seed, centr_region_pu_sub, et_3by3_cntr, et_jet_boosted, rEta_jet_boosted, rPhi_jet_boosted);

	// e-gamma algorithm
	egamma(algo_config.egamma_seed, algo_config.egamma_IsoFact, centr_region_pu_sub, et_3by3_cntr, nonIso_egamma_et, Iso_egamma_et);

	// Tau algorithm
	tau(algo_config.tau_seed, algo_config.tau_IsoFact, centr_region_pu_sub, et_3by3_cntr, nonIso_tau_et, Iso_tau_et);

////////////////////////////////////////////////////////////
	// Prepare algorithm results for sorting
	// so_in_x : sort object input
	// so_out_x: sort object output
	for (unsigned int phi = 0; phi < 18; phi++)
	{
#pragma HLS UNROLL
		for (int reg = 0; reg < 26; reg++)
		{
#pragma HLS UNROLL
			int idx_in = 26 * phi + reg;

			if (reg <= 6)
			{
				int idx_out = 14 * phi + reg;

				so_in_jet_fwd[idx_out].et = et_jet_calo[idx_in];
				so_in_jet_fwd[idx_out].idx = idx_in;
				so_in_jet_fwd[idx_out].rloc_phi = calo_regions[idx_in].rloc_phi;
				so_in_jet_fwd[idx_out].rloc_eta = calo_regions[idx_in].rloc_eta;
			}
			else if (reg >= 19)
			{
				int idx_out = 14 * phi - 12 + reg;

				so_in_jet_fwd[idx_out].et = et_jet_calo[idx_in];
				so_in_jet_fwd[idx_out].idx = idx_in;
				so_in_jet_fwd[idx_out].rloc_phi = calo_regions[idx_in].rloc_phi;
				so_in_jet_fwd[idx_out].rloc_eta = calo_regions[idx_in].rloc_eta;
			}
			else
			{
				int idx_out = 12 * phi - 7 + reg;

				so_in_jet_cr[idx_out].et = et_jet_calo[idx_in];
				so_in_jet_cr[idx_out].idx = idx_in;
				so_in_jet_cr[idx_out].rloc_phi = calo_regions[idx_in].rloc_phi;
				so_in_jet_cr[idx_out].rloc_eta = calo_regions[idx_in].rloc_eta;
			}
		}
	}

	for (int idx = 252; idx < 256; idx++)
	{
#pragma HLS UNROLL

		so_in_jet_fwd[idx].et = 0;
		so_in_jet_fwd[idx].idx = 0;
		so_in_jet_fwd[idx].rloc_phi = 0;
		so_in_jet_fwd[idx].rloc_eta = 0;
	}

	for (int idx = 216; idx < 256; idx++)
	{
#pragma HLS UNROLL

		so_in_jet_cr[idx].et = 0;
		so_in_jet_cr[idx].idx = 0;
		so_in_jet_cr[idx].rloc_phi = 0;
		so_in_jet_cr[idx].rloc_eta = 0;
	}

	for (int idx = 0; idx < 252; idx++)
	{
#pragma HLS UNROLL

		so_in_eg_noniso[idx].et = nonIso_egamma_et[idx];
		so_in_eg_noniso[idx].idx = idx;
		so_in_eg_noniso[idx].rloc_phi = centr_region[idx].rloc_phi;
		so_in_eg_noniso[idx].rloc_eta = centr_region[idx].rloc_eta;

		so_in_eg_iso[idx].et = Iso_egamma_et[idx];
		so_in_eg_iso[idx].idx = idx;
		so_in_eg_iso[idx].rloc_phi = centr_region[idx].rloc_phi;
		so_in_eg_iso[idx].rloc_eta = centr_region[idx].rloc_eta;

		so_in_tau_noniso[idx].et = nonIso_tau_et[idx];
		so_in_tau_noniso[idx].idx = idx;
		so_in_tau_noniso[idx].rloc_phi = centr_region[idx].rloc_phi;
		so_in_tau_noniso[idx].rloc_eta = centr_region[idx].rloc_eta;

		so_in_tau_iso[idx].et = Iso_tau_et[idx];
		so_in_tau_iso[idx].idx = idx;
		so_in_tau_iso[idx].rloc_phi = centr_region[idx].rloc_phi;
		so_in_tau_iso[idx].rloc_eta = centr_region[idx].rloc_eta;
                         
		so_in_jet_boosted[idx].et = et_jet_boosted[idx];
		so_in_jet_boosted[idx].idx = idx;
		so_in_jet_boosted[idx].rloc_phi = centr_region[idx].rloc_phi;
		so_in_jet_boosted[idx].rloc_eta = centr_region[idx].rloc_eta;

	}

	for (int idx = 252; idx < 256; idx++)
	{
#pragma HLS UNROLL

		so_in_eg_noniso[idx].et = 0;
		so_in_eg_noniso[idx].idx = 0;
		so_in_eg_noniso[idx].rloc_phi = 0;
		so_in_eg_noniso[idx].rloc_eta = 0;

		so_in_eg_iso[idx].et = 0;
		so_in_eg_iso[idx].idx = 0;
		so_in_eg_iso[idx].rloc_phi = 0;
		so_in_eg_iso[idx].rloc_eta = 0;

		so_in_tau_noniso[idx].et = 0;
		so_in_tau_noniso[idx].idx = 0;
		so_in_tau_noniso[idx].rloc_phi = 0;
		so_in_tau_noniso[idx].rloc_eta = 0;

		so_in_tau_iso[idx].et = 0;
		so_in_tau_iso[idx].idx = 0;
		so_in_tau_iso[idx].rloc_phi = 0;
		so_in_tau_iso[idx].rloc_eta = 0;

		so_in_jet_boosted[idx].et = 0;
		so_in_jet_boosted[idx].idx = 0;
		so_in_jet_boosted[idx].rloc_phi = 0;
		so_in_jet_boosted[idx].rloc_eta = 0;
	}

////////////////////////////////////////////////////////////
	// Do the actual sorting
	am_sort_256x8(so_in_jet_cr, so_out_jet_cr);
	am_sort_256x8(so_in_jet_fwd, so_out_jet_fwd);
	am_sort_256x8(so_in_jet_boosted, so_out_jet_boosted);

	am_sort_256x8(so_in_eg_noniso, so_out_eg_noniso);
	am_sort_256x8(so_in_eg_iso, so_out_eg_iso);

	am_sort_256x8(so_in_tau_noniso, so_out_tau_noniso);
	am_sort_256x8(so_in_tau_iso, so_out_tau_iso);

////////////////////////////////////////////////////////////
	// Assign the top 8 candidates to algorithnm outputs
	int olink;
	for (int idx = 0; idx < 8; idx++)
	{
#pragma HLS UNROLL
		ap_uint<1> side;
		ap_uint<9> idx_srt;
		ap_uint<6> ieta;


		{ // Central + Boosted Jets
			int bLo1 = 16;
			int bHi1 = bLo1 + 10; // 11 bits
			tmp_link_out[idx].range(bHi1, bLo1) = so_out_jet_cr[idx].et;
			idx_srt = so_out_jet_cr[idx].idx;

			side = calo_coor_full[idx_srt].side;
			int bLo2 = bHi1 + 1;
			int bHi2 = bLo2; // 1 bit
			tmp_link_out[idx].range(bHi2, bLo2) = side;

			int bLo3 = bHi2 + 1;
			int bHi3 = bLo3 + 6; // 7 bits
			tmp_link_out[idx].range(bHi3, bLo3) = calo_coor_full[idx_srt].iphi + so_out_jet_cr[idx].rloc_phi;

			int bLo4 = bHi3 + 1;
			int bHi4 = bLo4 + 5; // 6 bits
			if (side == 1)
				tmp_link_out[idx].range(bHi4, bLo4) = calo_coor_full[idx_srt].ieta - so_out_jet_cr[idx].rloc_eta;
			else
				tmp_link_out[idx].range(bHi4, bLo4) = calo_coor_full[idx_srt].ieta + so_out_jet_cr[idx].rloc_eta;

			if((double)tmp_link_out[idx].range(bHi1, bLo1) > 0) cout << "Jet CR " << idx << " ET: " << dec << tmp_link_out[idx].range(bHi1, bLo1) << " Side: " << tmp_link_out[idx].range(bHi2, bLo2) << " iPhi: " << tmp_link_out[idx].range(bHi3, bLo3)  << " iEta: " << tmp_link_out[idx].range(bHi4, bLo4) << endl ;
		}

		{ // Forward Jets
			int bLo5 = 16 + 25 + 1;
			int bHi5 = bLo5 + 10;
			tmp_link_out[idx].range(bHi5, bLo5) = so_out_jet_fwd[idx].et;
			idx_srt = so_out_jet_fwd[idx].idx;

			side = calo_coor_full[idx_srt].side;
			int bLo6 = bHi5 + 1;
			int bHi6 = bLo6;
			tmp_link_out[idx].range(bHi6, bLo6) = side;

			int bLo7 = bHi6 + 1;
			int bHi7 = bLo7 + 6; 
			tmp_link_out[idx].range(bHi7, bLo7) = calo_coor_full[idx_srt].iphi + so_out_jet_fwd[idx].rloc_phi;

			int bLo8 = bHi7 + 1;
			int bHi8 = bLo8 + 5;
			if (side == 1)
				tmp_link_out[idx].range(bHi8, bLo8) = calo_coor_full[idx_srt].ieta - so_out_jet_fwd[idx].rloc_eta;
			else
				tmp_link_out[idx].range(bHi8, bLo8) = calo_coor_full[idx_srt].ieta + so_out_jet_fwd[idx].rloc_eta;

			if((double)tmp_link_out[idx].range(bHi5, bLo5) > 0) cout << "Jet FWD " << idx << " ET: " << dec << tmp_link_out[idx].range(bHi5, bLo5) << " Side: " << tmp_link_out[idx].range(bHi6, bLo6) << " iPhi: " << tmp_link_out[idx].range(bHi7, bLo7)  << " iEta: " << tmp_link_out[idx].range(bHi8, bLo8) << endl ;
		}
		
		{ // Boosted jets
			int bLo9 = 16 + 25 + 25 + 1;
			int bHi9 = bLo9 + 10;
			tmp_link_out[idx].range(bHi9, bLo9) = so_out_jet_boosted[idx].et;
			idx_srt = so_out_jet_boosted[idx].idx;
			
			side = calo_coor[idx_srt].side;
			int bLo10 = bHi9 + 1;
			int bHi10 = bLo10;
			tmp_link_out[idx].range(bHi10, bLo10) = side;
	
			int bLo11 = bHi10 + 1;
			int bHi11 = bLo11 + 6;
			tmp_link_out[idx].range(bHi11, bLo11) = calo_coor[idx_srt].iphi + so_out_jet_boosted[idx].rloc_phi;

			int bLo12 = bHi11 + 1;
			int bHi12 = bLo12 + 5;
			if (side == 1)
				tmp_link_out[idx].range(bHi12, bLo12) = calo_coor[idx_srt].ieta - so_out_jet_boosted[idx].rloc_eta;
			else
				tmp_link_out[idx].range(bHi12, bLo12) = calo_coor[idx_srt].ieta + so_out_jet_boosted[idx].rloc_eta;

			int bLo13 = bHi12 + 1;
			int bHi13 = bLo13 + 2;
			tmp_link_out[idx].range(bHi13, bLo13) = rPhi_jet_boosted[idx_srt].to_ulong();

			int bLo14 = bHi13 + 1;
			int bHi14 = bLo14 + 2;
			tmp_link_out[idx].range(bHi14, bLo14) = rEta_jet_boosted[idx_srt].to_ulong();

			//if((double)tmp_link_out[idx].range(bHi9, bLo9) > 0) cout << "Jet Boosted " << idx << " ET: " << dec << tmp_link_out[idx].range(bHi9, bLo9) << " Side: " << tmp_link_out[idx].range(bHi10, bLo10) << " iPhi: " << tmp_link_out[idx].range(bHi11, bLo11)  << " iEta: " << tmp_link_out[idx].range(bHi12, bLo12) << " rPhi: " << tmp_link_out[idx].range(bHi13, bLo13).to_ulong() << " (" << rPhi_jet_boosted[idx_srt].to_string() << ") " << " rEta: " << tmp_link_out[idx].range(bHi14, bLo14).to_ulong() << " (" << rEta_jet_boosted[idx_srt].to_string() << ") " << endl ; 
			if((double)tmp_link_out[idx].range(bHi9, bLo9) > 0) cout << "Jet Boosted " << idx << " ET: " << dec << tmp_link_out[idx].range(bHi9, bLo9) << " Side: " << tmp_link_out[idx].range(bHi10, bLo10) << " iPhi: " << tmp_link_out[idx].range(bHi11, bLo11)  << " iEta: " << tmp_link_out[idx].range(bHi12, bLo12) << " rPhi: " << tmp_link_out[idx].range(bHi13, bLo13).to_ulong() << " rEta: " << tmp_link_out[idx].range(bHi14, bLo14).to_ulong() << endl;
		}

		{ // EG NonIso
			int bLo15 = 16 + 25 + 25 + 31 + 1;
			int bHi15 = bLo15 + 8;
			tmp_link_out[idx].range(bHi15, bLo15) = so_out_eg_noniso[idx].et;
			idx_srt = so_out_eg_noniso[idx].idx;

			side = calo_coor[idx_srt].side;
			int bLo16 = bHi15 + 1;
			int bHi16 = bLo16;
			tmp_link_out[idx].range(bHi16, bLo16) = side;

			int bLo17 = bHi16 + 1;
			int bHi17 = bLo17 + 6;
			tmp_link_out[idx].range(bHi17, bLo17) = calo_coor[idx_srt].iphi + so_out_eg_noniso[idx].rloc_phi;

			int bLo18 = bHi17 + 1;
			int bHi18 = bLo18 + 5;
			if (side == 1)
				tmp_link_out[idx].range(bHi18, bLo18) = calo_coor[idx_srt].ieta - so_out_eg_noniso[idx].rloc_eta;
			else
				tmp_link_out[idx].range(bHi18, bLo18) = calo_coor[idx_srt].ieta + so_out_eg_noniso[idx].rloc_eta;

			if((double)tmp_link_out[idx].range(bHi15, bLo15) > 0) cout << "EG NonIso " << idx << " ET: " << dec << tmp_link_out[idx].range(bHi15, bLo15) << " Side: " << tmp_link_out[idx].range(bHi16, bLo16) << " iPhi: " << tmp_link_out[idx].range(bHi17, bLo17)  << " iEta: " << tmp_link_out[idx].range(bHi18, bLo18) << endl ;
		}

		{ //EG ISO
			int bLo19 = 16 + 25 + 25 + 31 + 23 + 1;
			int bHi19 = bLo19 + 8;
			tmp_link_out[idx].range(bHi19, bLo19) = so_out_eg_iso[idx].et;
			idx_srt = so_out_eg_iso[idx].idx;

			side = calo_coor[idx_srt].side;
			int bLo20 = bHi19 + 1;
			int bHi20 = bLo20;
			tmp_link_out[idx].range(bHi20, bLo20) = side;

			int bLo21 = bHi20 + 1;
			int bHi21 = bLo21 + 6;
			tmp_link_out[idx].range(bHi21, bLo21) = calo_coor[idx_srt].iphi + so_out_eg_iso[idx].rloc_phi;

			int bLo22 = bHi21 + 1;
			int bHi22 = bLo22 + 5;
			if (side == 1)
				tmp_link_out[idx].range(bHi22, bLo22) = calo_coor[idx_srt].ieta - so_out_eg_iso[idx].rloc_eta;
			else
				tmp_link_out[idx].range(bHi22, bLo22) = calo_coor[idx_srt].ieta + so_out_eg_iso[idx].rloc_eta;

			if((double)tmp_link_out[idx].range(bHi19, bLo19) > 0) cout << "EG Iso " << idx << " ET: " << dec << tmp_link_out[idx].range(bHi19, bLo19) << " Side: " << tmp_link_out[idx].range(bHi20, bLo20) << " iPhi: " << tmp_link_out[idx].range(bHi21, bLo21)  << " iEta: " << tmp_link_out[idx].range(bHi22, bLo22) << endl ;
		}

		{ // Tau NonIso
			int bLo23 = 16 + 25 + 25 + 31 + 23 + 23 + 1;
			int bHi23 = bLo23 + 8;
			tmp_link_out[idx].range(bHi23, bLo23) = so_out_tau_noniso[idx].et;
			idx_srt = so_out_tau_noniso[idx].idx;

			side = calo_coor[idx_srt].side;
			int bLo24 = bHi23 + 1;
			int bHi24 = bLo24;
			tmp_link_out[idx].range(bHi24, bLo24) = side;

			int bLo25 = bHi24 + 1;
			int bHi25 = bLo25 + 6;
			tmp_link_out[idx].range(bHi25, bLo25) = calo_coor[idx_srt].iphi + so_out_tau_noniso[idx].rloc_phi;

			int bLo26 = bHi25 + 1;
			int bHi26 = bLo26 + 5;
			if (side == 1)
				tmp_link_out[idx].range(bHi26, bLo26) = calo_coor[idx_srt].ieta - so_out_tau_noniso[idx].rloc_eta;
			else
				tmp_link_out[idx].range(bHi26, bLo26) = calo_coor[idx_srt].ieta + so_out_tau_noniso[idx].rloc_eta;

			if((double)tmp_link_out[idx].range(bHi23, bLo23) > 0) cout << "Tau NonIso " << idx << " ET: " << dec << tmp_link_out[idx].range(bHi23, bLo23) << " Side: " << tmp_link_out[idx].range(bHi24, bLo24) << " iPhi: " << tmp_link_out[idx].range(bHi25, bLo25)  << " iEta: " << tmp_link_out[idx].range(bHi26, bLo26) << endl ;
		}

		{ //Tau ISO
			int bLo27 = 16 + 25 + 25 + 31 + 23 + 23 + 23 + 1;
			int bHi27 = bLo27 + 8;
			tmp_link_out[idx].range(bHi27, bLo27) = so_out_tau_iso[idx].et;
			idx_srt = so_out_tau_iso[idx].idx;

			side = calo_coor[idx_srt].side;
			int bLo28 = bHi27 + 1;
			int bHi28 = bLo28;
			tmp_link_out[idx].range(bHi28, bLo28) = side;

			int bLo29 = bHi28 + 1;
			int bHi29 = bLo29 + 6;
			tmp_link_out[idx].range(bHi29, bLo29) = calo_coor[idx_srt].iphi + so_out_tau_iso[idx].rloc_phi;

			int bLo30 = bHi29 + 1;
			int bHi30 = bLo30 + 5;
			if (side == 1)
				tmp_link_out[idx].range(bHi30, bLo30) = calo_coor[idx_srt].ieta - so_out_tau_iso[idx].rloc_eta;
			else
				tmp_link_out[idx].range(bHi30, bLo30) = calo_coor[idx_srt].ieta + so_out_tau_iso[idx].rloc_eta;

			if((double)tmp_link_out[idx].range(bHi27, bLo27) > 0) cout << "Tau Iso " << idx << " ET: " << dec << tmp_link_out[idx].range(bHi27, bLo27) << " Side: " << tmp_link_out[idx].range(bHi28, bLo28) << " iPhi: " << tmp_link_out[idx].range(bHi29, bLo29)  << " iEta: " << tmp_link_out[idx].range(bHi30, bLo30) << endl ;
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

