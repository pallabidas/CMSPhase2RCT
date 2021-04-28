#ifndef UCT_SUMMARY_CARD_H_
#define UCT_SUMMARY_CARD_H_

#include <stdint.h>
#include <ap_int.h>
#include <bitset>

using std::bitset;

// Number of central regions
#define NR_CNTR_REG (252)
// Number of forward regions
#define NR_FWD_REG (216)
// Number of calorimeter regions (central + forward)
#define NR_CALO_REG (NR_CNTR_REG + NR_FWD_REG)
// Number of central super-regions (2x2)
#define NR_SCNTR_REG (63)

#define PUM_LEVEL_BITSIZE (9)

typedef struct
{
	ap_uint<10> pum_thr;
	ap_uint<10> jet_seed;

	ap_uint<10> egamma_seed;
	ap_ufixed<7, 1, AP_RND, AP_SAT> egamma_IsoFact;

	ap_uint<10> tau_seed;
	ap_ufixed<7, 1, AP_RND, AP_SAT> tau_IsoFact;

} algo_config_t;

typedef struct
{
	ap_uint<10> et;
	bool eg_veto;
	bool tau_veto;
	ap_uint<2> rloc_phi;
	ap_uint<2> rloc_eta;
} region_t;

ap_uint<8> popcount(ap_uint<NR_CNTR_REG> bitString);

void et_3by3(ap_uint<10> et[NR_CNTR_REG], ap_uint<10> et_3by3[NR_CNTR_REG]);

void boostedjet(ap_uint<10> jet_seed, 
		region_t regions[NR_CNTR_REG], 
		ap_uint<10> et_3by3[NR_CNTR_REG], 
		ap_uint<10> jet_et[NR_SCNTR_REG], 
		bitset<3> rEta_jet[NR_SCNTR_REG], 
		bitset<3> rPhi_jet[NR_SCNTR_REG],
		ap_uint<9> rIdx_boostedjet[NR_SCNTR_REG]);

#endif
