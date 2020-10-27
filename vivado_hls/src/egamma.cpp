#include <cstdlib>
#include "ap_int.h"

#include "UCTSummaryCard.hpp"
#include "region_neighbors.h"

void egamma(ap_uint<10> egamma_seed,
		ap_ufixed<7, 1, AP_RND, AP_SAT> egamma_IsoFact,
		region_t regions[NR_CNTR_REG], ap_uint<10> et_3by3[NR_CNTR_REG],
		ap_uint<10> nonIso_egamma_et[NR_CNTR_REG],
		ap_uint<10> Iso_egamma_et[NR_CNTR_REG])
{
#pragma HLS inline
#pragma HLS PIPELINE II=3

#pragma HLS INTERFACE ap_none port=egamma_IsoFact
#pragma HLS INTERFACE ap_none port=egamma_seed

#pragma HLS ARRAY_RESHAPE variable=regions complete dim=1
#pragma HLS ARRAY_RESHAPE variable=et_3by3 complete dim=1
#pragma HLS ARRAY_RESHAPE variable=nonIso_egamma_et complete dim=1
#pragma HLS ARRAY_RESHAPE variable=Iso_egamma_et complete dim=1

	label0: for (int idx = 0; idx < NR_CNTR_REG; idx++)
	{
#pragma HLS UNROLL

		ap_uint<10> et_egamma;
		ap_uint<10> et_egamma_Iso;

		if ((idx % 14 == 0) || ((idx + 1) % 14) == 0)
		{
			et_egamma = 0;
			et_egamma_Iso = 0;
		}
		else
		{
			et_egamma = regions[idx].et;

			ap_int<10> neigh_N = region_cntr_neighbors[idx].nb_N;
			ap_int<10> neigh_S = region_cntr_neighbors[idx].nb_S;
			ap_int<10> neigh_W = region_cntr_neighbors[idx].nb_W;
			ap_int<10> neigh_E = region_cntr_neighbors[idx].nb_E;

			if (regions[idx].eg_veto == false && regions[idx].et > egamma_seed)
			{
/////// North
				if (neigh_N != -1)
				{
					if (regions[idx].rloc_phi == 2
							&& regions[neigh_N].rloc_phi == 0
							&& std::abs(
									(ap_int<3> ) regions[idx].rloc_eta
											- (ap_int<3> ) regions[neigh_N].rloc_eta)
									< 2 && regions[neigh_N].eg_veto == false)
					{
						if (regions[idx].et >= regions[neigh_N].et)
						{
							et_egamma += regions[neigh_N].et;
						}
						else
						{
							et_egamma = 0;
						}
					}
				}
/////// South
				if (neigh_S != -1)
				{
					if (regions[idx].rloc_phi == 0
							&& regions[neigh_S].rloc_phi == 2
							&& std::abs(
									(ap_int<3> ) regions[idx].rloc_eta
											- (ap_int<3> ) regions[neigh_S].rloc_eta)
									< 2 && regions[neigh_S].eg_veto == false)
					{
						if (regions[idx].et > regions[neigh_S].et)
						{
							et_egamma += regions[neigh_S].et;
						}
						else
						{
							et_egamma = 0;
						}
					}
				}

/////// West
				if (neigh_W != -1)
				{
					if (regions[idx].rloc_phi == 0
							&& regions[neigh_W].rloc_phi == 2
							&& std::abs(
									(ap_int<3> ) regions[idx].rloc_phi
											- (ap_int<3> ) regions[neigh_W].rloc_eta)
									< 2 && regions[neigh_W].eg_veto == false)
					{
						if (regions[idx].et >= regions[neigh_W].et)
						{
							et_egamma += regions[neigh_W].et;
						}
						else
						{
							et_egamma = 0;
						}
					}
				}
				/////// East
				if (neigh_E != -1)
				{
					if (regions[idx].rloc_phi == 2
							&& regions[neigh_E].rloc_phi == 0
							&& std::abs(
									(ap_int<3> ) regions[idx].rloc_phi
											- (ap_int<3> ) regions[neigh_E].rloc_eta)
									< 2 && regions[neigh_E].eg_veto == false)
					{
						if (regions[idx].et > regions[neigh_E].et)
						{
							et_egamma += regions[neigh_E].et;
						}
						else
						{
							et_egamma = 0;
						}
					}
				}

			}
			else
			{
				et_egamma = 0;
			}
			ap_uint<10> isolation;

			if (et_3by3[idx] > et_egamma)
				isolation = et_3by3[idx] - et_egamma;
			else
				isolation = 0;

			if (isolation < (ap_uint<10> ) (egamma_IsoFact * et_egamma))
				et_egamma_Iso = et_egamma;
			else
				et_egamma_Iso = 0;
		}

		nonIso_egamma_et[idx] = et_egamma;
		Iso_egamma_et[idx] = et_egamma_Iso;

	}
	return;
}
