#ifndef __AM_SORT_H__
#define __AM_SORT_H__

#include "ap_int.h"

// Sort Object (SO)
typedef struct
{
	ap_uint<10> et;
	ap_uint<9> idx;
	ap_uint<2> rloc_phi;
	ap_uint<2> rloc_eta;
} t_so;

void bm2_plus_mod(t_so so_in[2], t_so & so_out1, t_so & so_out2 );

void bm2_plus(t_so so_in[2], t_so so_out[2]);
void bm2_minus(t_so so_in[2], t_so so_out[2]);
void bm4_plus(t_so so_in[4], t_so so_out[4]);
void bm4_minus(t_so so_in[4], t_so so_out[4]);
void bm8_plus(t_so so_in[8], t_so so_out[8]);
void bm16_8_plus(t_so so_in[16], t_so so_out[8]);

void oem4(t_so so_in[4], t_so so_out[4]);
void oem8(t_so so_in[8], t_so so_out[8]);
void oem16(t_so so_in[16], t_so so_out[16]);
void oem16_8(t_so so_in[16], t_so so_out[8]);

void max16(t_so so_in[16], t_so so_out[8]);

void am_sort_256x8(t_so so_in[256], t_so so_out[8]);

void oem8_f(t_so so_in[8], t_so so_out[8]);

#endif
