#include "UCTSummaryCard.hpp"

#include "adder_tree.h"
#include "am_sort.h"

#include "PU_LUT.h"
#include "calo_out_coordinates.h"

void UCTSummaryCard(
		region_t centr_region[NR_CNTR_REG],  // central region inputs
		region_t fwd_region[NR_FWD_REG],     // forward region inputs
		algo_config_t algo_config,           // algoritm configuration
		algo_outputs_t & algo_outputs        // algorithm outputs
		)
{

	region_t calo_regions[NR_CALO_REG];

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

	region_t centr_region_pu_sub[NR_CNTR_REG];

////////////////////////////////////////////////////////////////////////

#pragma HLS INTERFACE ap_none port=algo_config

#pragma HLS PIPELINE II=6

#pragma HLS ARRAY_RESHAPE variable=calo_regions complete dim=1
#pragma HLS ARRAY_RESHAPE variable=centr_region complete dim=1
#pragma HLS ARRAY_RESHAPE variable=fwd_region complete dim=1
#pragma HLS ARRAY_RESHAPE variable=centr_region_pu_sub complete dim=1

#pragma HLS ARRAY_RESHAPE variable=algo_outputs.jet_central complete dim=1
#pragma HLS ARRAY_RESHAPE variable=algo_outputs.jet_forward complete dim=1
#pragma HLS ARRAY_RESHAPE variable=algo_outputs.jet_boosted complete dim=1
#pragma HLS ARRAY_RESHAPE variable=algo_outputs.eg_noniso complete dim=1
#pragma HLS ARRAY_RESHAPE variable=algo_outputs.eg_iso complete dim=1
#pragma HLS ARRAY_RESHAPE variable=algo_outputs.tau_noniso complete dim=1
#pragma HLS ARRAY_RESHAPE variable=algo_outputs.tau_iso complete dim=1

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
				calo_regions[26 * phi + reg] = fwd_region[12 * phi + reg];
			}
			else if (reg >= 20) // FWD+
			{
				calo_regions[26 * phi + reg] = fwd_region[12 * phi - 14 + reg];
			}
			else // CNTR- and CNTR+
			{
				calo_regions[26 * phi + reg] = centr_region[14 * phi - 6 + reg];
			}
		}
	}

////////////////////////////////////////////////////////////
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

	algo_outputs.et.et = sums.et_sum;
	algo_outputs.ht.et = sums.ht_sum;
	algo_outputs.met.et = sums.met_sum;
	algo_outputs.mht.et = sums.ht_sum;

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
	//boostedjet(algo_config.jet_seed, calo_regions, pu_sub_et_calo, et_3by3_calo, et_jet_boosted, rEta_jet_boosted, rPhi_jet_boosted);
	boostedjet(algo_config.jet_seed, centr_region_pu_sub, et_3by3_cntr, et_jet_boosted, rEta_jet_boosted, rPhi_jet_boosted);

	// e-gamma algorithm
	egamma(algo_config.egamma_seed, algo_config.egamma_IsoFact,
			centr_region_pu_sub, et_3by3_cntr, nonIso_egamma_et, Iso_egamma_et);

	// Tau algorithm
	tau(algo_config.tau_seed, algo_config.tau_IsoFact, centr_region_pu_sub,
			et_3by3_cntr, nonIso_tau_et, Iso_tau_et);

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

				//so_in_jet_boosted[idx_out].et = et_jet_boosted[idx_in];
				//so_in_jet_boosted[idx_out].idx = idx_in;
				//so_in_jet_boosted[idx_out].rloc_phi = rEta_jet_boosted[idx_in];
				//so_in_jet_boosted[idx_out].rloc_eta = rPhi_jet_boosted[idx_in];
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

		//so_in_jet_boosted[idx].et = 0;
		//so_in_jet_boosted[idx].idx = 0;
		//so_in_jet_boosted[idx].rloc_phi = 0;
		//so_in_jet_boosted[idx].rloc_eta = 0;
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
	for (int idx = 0; idx < 8; idx++)
	{
#pragma HLS UNROLL
		ap_uint<1> side;
		ap_uint<9> idx_srt;
		ap_uint<6> ieta;

		{ // Central Jets
			algo_outputs.jet_central[idx].et = so_out_jet_cr[idx].et;
			idx_srt = so_out_jet_cr[idx].idx;

			side = calo_coor_full[idx_srt].side;
			algo_outputs.jet_central[idx].side = side;

			algo_outputs.jet_central[idx].iphi = calo_coor_full[idx_srt].iphi
					+ so_out_jet_cr[idx].rloc_phi;

			if (side == 1)
				algo_outputs.jet_central[idx].ieta =
						calo_coor_full[idx_srt].ieta
								- so_out_jet_cr[idx].rloc_eta;
			else
				algo_outputs.jet_central[idx].ieta =
						calo_coor_full[idx_srt].ieta
								+ so_out_jet_cr[idx].rloc_eta;
		}

		{ // Forward Jets
			algo_outputs.jet_forward[idx].et = so_out_jet_fwd[idx].et;
			idx_srt = so_out_jet_fwd[idx].idx;

			side = calo_coor_full[idx_srt].side;
			algo_outputs.jet_forward[idx].side = side;

			algo_outputs.jet_forward[idx].iphi = calo_coor_full[idx_srt].iphi
					+ so_out_jet_fwd[idx].rloc_phi;

			if (side == 1)
				algo_outputs.jet_forward[idx].ieta =
						calo_coor_full[idx_srt].ieta
								- so_out_jet_fwd[idx].rloc_eta;
			else
				algo_outputs.jet_forward[idx].ieta =
						calo_coor_full[idx_srt].ieta
								+ so_out_jet_fwd[idx].rloc_eta;
		}
		
		{ // Boosted jets
			algo_outputs.jet_boosted[idx].et = so_out_jet_boosted[idx].et;
			idx_srt = so_out_jet_boosted[idx].idx;
			
			side = calo_coor[idx_srt].side;
			algo_outputs.jet_boosted[idx].side = side;
	
			algo_outputs.jet_boosted[idx].iphi = calo_coor[idx_srt].iphi
					+ so_out_jet_boosted[idx].rloc_phi;

			if (side == 1)
				algo_outputs.jet_boosted[idx].ieta =
						calo_coor[idx_srt].ieta
								- so_out_jet_boosted[idx].rloc_eta;
			else
				algo_outputs.jet_boosted[idx].ieta =
						calo_coor[idx_srt].ieta
								+ so_out_jet_boosted[idx].rloc_eta;

			algo_outputs.jet_boosted[idx].rEta = rEta_jet_boosted[idx_srt];
			algo_outputs.jet_boosted[idx].rPhi = rPhi_jet_boosted[idx_srt];

		}

		{ // EG NonIso
			algo_outputs.eg_noniso[idx].et = so_out_eg_noniso[idx].et;
			idx_srt = so_out_eg_noniso[idx].idx;

			side = calo_coor[idx_srt].side;
			algo_outputs.eg_noniso[idx].side = side;

			algo_outputs.eg_noniso[idx].iphi = calo_coor[idx_srt].iphi
					+ so_out_eg_noniso[idx].rloc_phi;

			if (side == 1)
				algo_outputs.eg_noniso[idx].ieta = calo_coor[idx_srt].ieta
						- so_out_eg_noniso[idx].rloc_eta;
			else
				algo_outputs.eg_noniso[idx].ieta = calo_coor[idx_srt].ieta
						+ so_out_eg_noniso[idx].rloc_eta;
		}

		{ //EG ISO
			algo_outputs.eg_iso[idx].et = so_out_eg_iso[idx].et;
			idx_srt = so_out_eg_iso[idx].idx;

			side = calo_coor[idx_srt].side;
			algo_outputs.eg_iso[idx].side = side;

			algo_outputs.eg_iso[idx].iphi = calo_coor[idx_srt].iphi
					+ so_out_eg_iso[idx].rloc_phi;

			if (side == 1)
				algo_outputs.eg_iso[idx].ieta = calo_coor[idx_srt].ieta
						- so_out_eg_iso[idx].rloc_eta;
			else
				algo_outputs.eg_iso[idx].ieta = calo_coor[idx_srt].ieta
						+ so_out_eg_iso[idx].rloc_eta;
		}

		{ // Tau NonIso
			algo_outputs.tau_noniso[idx].et = so_out_tau_noniso[idx].et;
			idx_srt = so_out_tau_noniso[idx].idx;

			side = calo_coor[idx_srt].side;
			algo_outputs.tau_noniso[idx].side = side;

			algo_outputs.tau_noniso[idx].iphi = calo_coor[idx_srt].iphi
					+ so_out_tau_noniso[idx].rloc_phi;

			if (side == 1)
				algo_outputs.tau_noniso[idx].ieta = calo_coor[idx_srt].ieta
						- so_out_tau_noniso[idx].rloc_eta;
			else
				algo_outputs.tau_noniso[idx].ieta = calo_coor[idx_srt].ieta
						+ so_out_tau_noniso[idx].rloc_eta;
		}

		{ //Tau ISO
			algo_outputs.tau_iso[idx].et = so_out_tau_iso[idx].et;
			idx_srt = so_out_tau_iso[idx].idx;

			side = calo_coor[idx_srt].side;
			algo_outputs.tau_iso[idx].side = side;

			algo_outputs.tau_iso[idx].iphi = calo_coor[idx_srt].iphi
					+ so_out_tau_iso[idx].rloc_phi;

			if (side == 1)
				algo_outputs.tau_iso[idx].ieta = calo_coor[idx_srt].ieta
						- so_out_tau_iso[idx].rloc_eta;
			else
				algo_outputs.tau_iso[idx].ieta = calo_coor[idx_srt].ieta
						+ so_out_tau_iso[idx].rloc_eta;
		}
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
