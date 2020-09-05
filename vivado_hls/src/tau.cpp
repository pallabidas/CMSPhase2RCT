#include <cstdlib>
#include "ap_int.h"

#include "UCTSummaryCard.hpp"
#include "region_neighbors.h"

void tau(ap_uint<10> tau_seed, ap_ufixed<7, 1, AP_RND, AP_SAT> tau_IsoFact,
		region_t regions[NR_CNTR_REG], ap_uint<10> et_3by3[NR_CNTR_REG],
		ap_uint<10> nonIso_tau_et[NR_CNTR_REG],
		ap_uint<10> Iso_tau_et[NR_CNTR_REG])
{
#pragma HLS PIPELINE II=3

#pragma HLS INTERFACE ap_none port=tau_IsoFact
#pragma HLS INTERFACE ap_none port=tau_seed

#pragma HLS ARRAY_RESHAPE variable=regions complete dim=1
#pragma HLS resource variable=regions core=AddSub_DSP
#pragma HLS ARRAY_RESHAPE variable=et_3by3 complete dim=1
#pragma HLS ARRAY_RESHAPE variable=nonIso_tau_et complete dim=1
#pragma HLS ARRAY_RESHAPE variable=Iso_tau_et complete dim=1

	label0: for (int idx = 0; idx < NR_CNTR_REG; idx++)
	{
#pragma HLS UNROLL
		ap_uint<10> et_tau;
		ap_uint<10> et_tau_Iso;

		if ((idx % 14 == 0) || ((idx + 1) % 14) == 0)
		{
			et_tau = 0;
			et_tau_Iso = 0;
		}
		else
		{
			//ap_uint<10> et_central = regions[idx].et;
			ap_uint<10> et_N = 0;
			ap_uint<10> et_S = 0;
			ap_uint<10> et_E = 0;
			ap_uint<10> et_W = 0;

			et_tau = regions[idx].et;

			ap_int<10> neigh_N = region_cntr_neighbors[idx].nb_N;
			ap_int<10> neigh_S = region_cntr_neighbors[idx].nb_S;
			ap_int<10> neigh_W = region_cntr_neighbors[idx].nb_W;
			ap_int<10> neigh_E = region_cntr_neighbors[idx].nb_E;

			if (regions[idx].tau_veto == false && regions[idx].et > tau_seed)
			{
				/////// North
				if (neigh_N != -1)
				{
					if (regions[idx].rloc_phi == 2
							&& regions[neigh_N].rloc_phi == 0
							&& std::abs(
									(ap_int<3> ) regions[idx].rloc_eta
											- (ap_int<3> ) regions[neigh_N].rloc_eta)
									< 2 && regions[neigh_N].tau_veto == false)
					{
						if (regions[idx].et >= regions[neigh_N].et)
						{
							et_tau += regions[neigh_N].et;
						}
						else
						{
							et_tau = 0;
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
									< 2 && regions[neigh_S].tau_veto == false)
					{
						if (regions[idx].et > regions[neigh_S].et)
						{
							et_tau += regions[neigh_S].et;
						}
						else
						{
							et_tau = 0;
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
									< 2 && regions[neigh_W].tau_veto == false)
					{
						if (regions[idx].et >= regions[neigh_W].et)
						{
							et_tau += regions[neigh_W].et;
						}
						else
						{
							et_tau = 0;
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
									< 2 && regions[neigh_E].tau_veto == false)
					{
						if (regions[idx].et > regions[neigh_E].et)
						{
							et_tau += regions[neigh_E].et;
						}
						else
						{
							et_tau = 0;
						}
					}
				}

			}
			else
			{
				et_tau = 0;
			}

			ap_uint<10> isolation;

			if (et_3by3[idx] > et_tau)
				isolation = et_3by3[idx] - et_tau;
			else
				isolation = 0;

			if (isolation < (ap_uint<10> ) (tau_IsoFact * et_tau))
				et_tau_Iso = et_tau;
			else
				et_tau_Iso = 0;
		}

		nonIso_tau_et[idx] = et_tau;
		Iso_tau_et[idx] = et_tau_Iso;

	}
	return;
}
