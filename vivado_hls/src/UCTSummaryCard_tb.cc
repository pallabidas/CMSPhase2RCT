
#include <stdio.h>
#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <bitset>
#include "UCTSummaryCard.hpp"

using namespace std;
using std::bitset;

vector< vector<unsigned int> >  read_regions_from_file(void)
		{
	ifstream inFile;
	inFile.open("data.txt");//open the input file

	stringstream strStream;
	strStream << inFile.rdbuf();//read the file
	string str = strStream.str();//str holds the content of the file

	string buf; // Have a buffer string
	stringstream ss(str); // Insert the string into a stream

	vector<string> tokens; // Create vector to hold our words

	while (ss >> buf)
		tokens.push_back(buf);

	unsigned int tokens_per_event = 2 + 18 * 26;
	unsigned int events_in_file = tokens.size() / tokens_per_event;

	vector< vector<unsigned int> > v;

	for (int i = 0; i < events_in_file; i++)
	{
		if (tokens[i * tokens_per_event] == "#Event")
		{
		//	cout << "Processing event #: " << tokens[i * tokens_per_event + 1] << endl;
		}
		vector<unsigned int> region;
		for (int j = 0; j < 18 * 26; j++)
		{
			unsigned int idx = i * tokens_per_event + 2 + j;

			std::stringstream ss2;
			unsigned int x;

			ss2 << hex << tokens[idx].c_str() ;
			ss2 >> x;

			region.push_back(x);
		}
		v.push_back(region);
	}
	return v;
		}

// FWD- regions        Central-           Central+              FWD+ regions
//0 1 2 3 4 5      6 7 8 9 10 11 12 13 14 15 16 17 18 19       20 21 22 23 24 25

void unpack_input_regions(vector<unsigned int>  region_raw,
		region_t centr_region[18 * 2 * 7],
		region_t fwd_region[18 * 2 * 6]
)
{
	for (unsigned int phi = 0; phi < 18; phi++)
	{
		for (int reg = 0; reg < 26; reg++)
		{
			if (reg <= 5)
			{
				unsigned int idx_in = 468 - ( phi + 1 ) * 26 + reg;
				unsigned int idx_out = phi * 12 + reg;

				fwd_region[idx_out].et = (region_raw[idx_in] & 0x3FF) >> 0;
				fwd_region[idx_out].eg_veto = (region_raw[idx_in] & 0x7FF) >> 10;
				fwd_region[idx_out].tau_veto = (region_raw[idx_in] & 0xFFF) >> 11;
				fwd_region[idx_out].rloc_phi = (region_raw[idx_in] & 0x3FFF) >> 12;
				fwd_region[idx_out].rloc_eta = (region_raw[idx_in] & 0xFFFF) >> 14;
			}
			else if (reg >= 20)
			{
				unsigned int idx_in = 468 - ( phi + 1 ) * 26 + reg;
				unsigned int idx_out = phi * 12 - 14 + reg;

				fwd_region[idx_out].et = (region_raw[idx_in] & 0x3FF) >> 0;
				fwd_region[idx_out].eg_veto = (region_raw[idx_in] & 0x7FF) >> 10;
				fwd_region[idx_out].tau_veto = (region_raw[idx_in] & 0xFFF) >> 11;
				fwd_region[idx_out].rloc_phi = (region_raw[idx_in] & 0x3FFF) >> 12;
				fwd_region[idx_out].rloc_eta = (region_raw[idx_in] & 0xFFFF) >> 14;

			}
			//if (reg > 5 && reg < 20)
			else

			{
				unsigned int idx_in = 468 - ( phi + 1 ) * 26 + reg;
				unsigned int idx_out = phi * 14 - 6 + reg;

				centr_region[idx_out].et =  (region_raw[idx_in] & 0x3FF >> 0);
				centr_region[idx_out].eg_veto = (region_raw[idx_in] & 0x7FF) >> 10;
				centr_region[idx_out].tau_veto = (region_raw[idx_in] & 0xFFF) >> 11;
				centr_region[idx_out].rloc_phi = (region_raw[idx_in] & 0x3FFF) >> 12;
				centr_region[idx_out].rloc_eta = (region_raw[idx_in] & 0xFFFF) >> 14;
			}
		}
	}

	return;
}

int main()
{

	region_t centr_region[NR_CNTR_REG]; // 252
	region_t fwd_region[NR_FWD_REG];   // 216

	cout << " ------ Start new event -----" << endl;

	vector< vector< unsigned int> > regions_raw = read_regions_from_file();

	unsigned int events_in_file = regions_raw.size();

	algo_outputs_t algo_outputs;
	algo_config_t algo_config;

	algo_config.egamma_IsoFact =  0.3;
	algo_config.egamma_seed = 5;
	algo_config.jet_seed = 10;
	algo_config.pum_thr = 0;
	algo_config.tau_IsoFact =  0.3;
	algo_config.tau_seed = 10;


	for (int i = 0; i < events_in_file; i++)
	{
		unpack_input_regions(regions_raw[i], centr_region, fwd_region );

		UCTSummaryCard(centr_region, fwd_region, algo_config, algo_outputs);


		cout << endl;
		cout << "----------------" << endl << "Event #: " << i << endl;

		for (int i = 0; i < NR_CNTR_REG; i++)
		{
			if (centr_region[i].et != 0)
			     cout << "Non-zero region.  idx: " << dec << i << " et: " << hex << centr_region[i].et << " Tau Veto: " <<  centr_region[i].tau_veto << " EG Veto: " <<  centr_region[i].eg_veto << " Phi: " <<  centr_region[i].rloc_phi << " Eta: " <<  centr_region[i].rloc_eta << endl;
		}

		cout << endl;
		cout << "----------------" << endl << "Outputs: " << endl << endl;


		for (int i = 0; i < 8; i++)
		{
			ap_uint<11> et = algo_outputs.jet_central[i].et;
			ap_uint<1> side = algo_outputs.jet_central[i].side;
			ap_uint<7> iphi = algo_outputs.jet_central[i].iphi;
			ap_uint<6> ieta = algo_outputs.jet_central[i].ieta;

			cout << "Jet CR" << i << " ET: " << hex << et << dec << " Side: " << side << " iPhi: " << iphi << " iEta: " << ieta << endl ;

		}

		cout << endl;

		for (int i = 0; i < 8; i++)
		{
			ap_uint<11> et = algo_outputs.jet_forward[i].et;
			ap_uint<1> side = algo_outputs.jet_forward[i].side;
			ap_uint<7> iphi = algo_outputs.jet_forward[i].iphi;
			ap_uint<6> ieta = algo_outputs.jet_forward[i].ieta;

			cout << "Jet FWD " << i << " ET: " << hex << et << dec << " Side: " << side << " iPhi: " << iphi << " iEta: " << ieta << endl ;

		}

		cout << endl;

		for (int i = 0; i < 8; i++)
		{
			ap_uint<11> et = algo_outputs.jet_boosted[i].et;
			ap_uint<1> side = algo_outputs.jet_boosted[i].side;
			ap_uint<7> iphi = algo_outputs.jet_boosted[i].iphi;
			ap_uint<6> ieta = algo_outputs.jet_boosted[i].ieta;
			bitset<3> rEta = algo_outputs.jet_boosted[i].rEta;
			bitset<3> rPhi = algo_outputs.jet_boosted[i].rPhi;

			cout << "Jet Boosted " << i << " ET: " << hex << et << dec << " Side: " << side << " iPhi: " << iphi << " iEta: " << ieta << " rEta: " << rEta.to_string() << " rPhi: " << rPhi.to_string() << endl;

		}

		cout << endl;

		for (int i = 0; i < 8; i++)
		{
			ap_uint<11> et = algo_outputs.eg_noniso[i].et;
			ap_uint<1> side = algo_outputs.eg_noniso[i].side;
			ap_uint<7> iphi = algo_outputs.eg_noniso[i].iphi;
			ap_uint<6> ieta = algo_outputs.eg_noniso[i].ieta;

			cout << "eg_noniso " << i << " ET: " << hex << et << dec << " Side: " << side << " iPhi: " << iphi << " iEta: " << ieta << endl ;

		}

		cout << endl;

		for (int i = 0; i < 8; i++)
		{
			ap_uint<11> et = algo_outputs.eg_iso[i].et;
			ap_uint<1> side = algo_outputs.eg_iso[i].side;
			ap_uint<7> iphi = algo_outputs.eg_iso[i].iphi;
			ap_uint<6> ieta = algo_outputs.eg_iso[i].ieta;

			cout << "eg_iso " << i << " ET: " << hex << et << dec << " Side: " << side << " iPhi: " << iphi << " iEta: " << ieta << endl ;

		}

		cout << endl;

		for (int i = 0; i < 8; i++)
		{
			ap_uint<11> et = algo_outputs.tau_noniso[i].et;
			ap_uint<1> side = algo_outputs.tau_noniso[i].side;
			ap_uint<7> iphi = algo_outputs.tau_noniso[i].iphi;
			ap_uint<6> ieta = algo_outputs.tau_noniso[i].ieta;

			cout << "tau_noniso " << i << " ET: " << hex << et << dec << " Side: " << side << " iPhi: " << iphi << " iEta: " << ieta << endl ;

		}

		cout << endl;

		for (int i = 0; i < 8; i++)
		{
			ap_uint<11> et = algo_outputs.tau_iso[i].et;
			ap_uint<1> side = algo_outputs.tau_iso[i].side;
			ap_uint<7> iphi = algo_outputs.tau_iso[i].iphi;
			ap_uint<6> ieta = algo_outputs.tau_iso[i].ieta;

			cout << "tau_iso " << i << " ET: " << hex << et << dec << " Side: " << side << " iPhi: " << iphi << " iEta: " << ieta << endl ;

		}

		cout << " ------ End new event -----" << endl;

	}

	return 0;
}

