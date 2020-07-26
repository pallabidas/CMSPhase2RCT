#include "ap_int.h"
#include "am_sort.h"


// Sorter network implementation, conceived by a former UW student Amin
// The implementation below was translated into HLS from his HDL source code

void am_sort_256x8(t_so so_in[256], t_so so_out[8])
{

#pragma HLS PIPELINE II=6

#pragma HLS ARRAY_RESHAPE variable=so_in complete dim=0
#pragma HLS ARRAY_RESHAPE variable=so_out complete dim=0

	t_so lvl0[128][2];
	t_so lvl1_2[128][2];
	t_so lvl1_4[64][4];
	t_so lvl3_4[64][4];
	t_so lvl3_8[32][8];
	t_so lvl6_8[32][8];
	t_so lvl6_16[16][16];
	t_so lvl10_8[16][8];
	t_so lvl10_16[8][16];
	t_so lvl14_8[8][8];
	t_so lvl14_16[4][16];
	t_so lvl18_8[4][8];
	t_so lvl18_16[2][16];
	t_so lvl22_8[2][8];
	t_so lvl22[16];
	t_so lvl23[8];

#pragma HLS ARRAY_RESHAPE variable=lvl0 complete dim=0
#pragma HLS ARRAY_RESHAPE variable=lvl1_2 complete dim=0
#pragma HLS ARRAY_RESHAPE variable=lvl1_4 complete dim=0
#pragma HLS ARRAY_RESHAPE variable=lvl3_4 complete dim=0
#pragma HLS ARRAY_RESHAPE variable=lvl3_8 complete dim=0
#pragma HLS ARRAY_RESHAPE variable=lvl6_8 complete dim=0
#pragma HLS ARRAY_RESHAPE variable=lvl6_16 complete dim=0
#pragma HLS ARRAY_RESHAPE variable=lvl10_8 complete dim=0
#pragma HLS ARRAY_RESHAPE variable=lvl10_16 complete dim=0
#pragma HLS ARRAY_RESHAPE variable=lvl14_8 complete dim=0
#pragma HLS ARRAY_RESHAPE variable=lvl14_16 complete dim=0
#pragma HLS ARRAY_RESHAPE variable=lvl18_8 complete dim=0
#pragma HLS ARRAY_RESHAPE variable=lvl18_16 complete dim=0
#pragma HLS ARRAY_RESHAPE variable=lvl22_8 complete dim=0
#pragma HLS ARRAY_RESHAPE variable=lvl22 complete dim=0
#pragma HLS ARRAY_RESHAPE variable=lvl23 complete dim=0

	for (int i = 0; i < 128; i++)
	{

//		lvl0[i][0] = so_in[i * 2];
//		lvl0[i][1] = so_in[i * 2 + 1];

		t_so tmp_in[2];
		t_so tmp_out[2];

#pragma HLS ARRAY_RESHAPE variable=tmp_in complete dim=0
#pragma HLS ARRAY_RESHAPE variable=tmp_out complete dim=0

		tmp_in[0] = so_in[i * 2];
		tmp_in[1] = so_in[i * 2 + 1];

		bm2_plus(tmp_in, tmp_out);

		lvl1_2[i][0] = tmp_out[0];
		lvl1_2[i][1] = tmp_out[1];

	}

	for (int i = 0; i < 64; i++)
	{

		t_so tmp_in[4];
		t_so tmp_out[4];

#pragma HLS ARRAY_RESHAPE variable=tmp_in complete dim=0
#pragma HLS ARRAY_RESHAPE variable=tmp_out complete dim=0

		tmp_in[0] = lvl1_2[i * 2][0];
		tmp_in[1] = lvl1_2[i * 2][1];
		tmp_in[2] = lvl1_2[i * 2 + 1][0];
		tmp_in[3] = lvl1_2[i * 2 + 1][1];

		oem4(tmp_in, tmp_out);

		lvl3_4[i][0] = tmp_out[0];
		lvl3_4[i][1] = tmp_out[1];
		lvl3_4[i][2] = tmp_out[2];
		lvl3_4[i][3] = tmp_out[3];

	}

	for (int i = 0; i < 32; i++)
	{

		t_so tmp_in[8];
		t_so tmp_out[8];

#pragma HLS ARRAY_RESHAPE variable=tmp_in complete dim=0
#pragma HLS ARRAY_RESHAPE variable=tmp_out complete dim=0

		tmp_in[0] = lvl3_4[i * 2][0];
		tmp_in[1] = lvl3_4[i * 2][1];
		tmp_in[2] = lvl3_4[i * 2][2];
		tmp_in[3] = lvl3_4[i * 2][3];
		tmp_in[4] = lvl3_4[i * 2 + 1][0];
		tmp_in[5] = lvl3_4[i * 2 + 1][1];
		tmp_in[6] = lvl3_4[i * 2 + 1][2];
		tmp_in[7] = lvl3_4[i * 2 + 1][3];

		oem8(tmp_in, tmp_out);

		lvl6_8[i][0] = tmp_out[0];
		lvl6_8[i][1] = tmp_out[1];
		lvl6_8[i][2] = tmp_out[2];
		lvl6_8[i][3] = tmp_out[3];
		lvl6_8[i][4] = tmp_out[4];
		lvl6_8[i][5] = tmp_out[5];
		lvl6_8[i][6] = tmp_out[6];
		lvl6_8[i][7] = tmp_out[7];

	}
	///

	for (int i = 0; i < 16; i++)
	{

		t_so tmp_in[16];
		t_so tmp_out[8];

#pragma HLS ARRAY_RESHAPE variable=tmp_in complete dim=0
#pragma HLS ARRAY_RESHAPE variable=tmp_out complete dim=0

		tmp_in[0] = lvl6_8[i * 2][0];
		tmp_in[1] = lvl6_8[i * 2][1];
		tmp_in[2] = lvl6_8[i * 2][2];
		tmp_in[3] = lvl6_8[i * 2][3];
		tmp_in[4] = lvl6_8[i * 2][4];
		tmp_in[5] = lvl6_8[i * 2][5];
		tmp_in[6] = lvl6_8[i * 2][6];
		tmp_in[7] = lvl6_8[i * 2][7];
		tmp_in[8] = lvl6_8[i * 2 + 1][0];
		tmp_in[9] = lvl6_8[i * 2 + 1][1];
		tmp_in[10] = lvl6_8[i * 2 + 1][2];
		tmp_in[11] = lvl6_8[i * 2 + 1][3];
		tmp_in[12] = lvl6_8[i * 2 + 1][4];
		tmp_in[13] = lvl6_8[i * 2 + 1][5];
		tmp_in[14] = lvl6_8[i * 2 + 1][6];
		tmp_in[15] = lvl6_8[i * 2 + 1][7];

		oem16_8(tmp_in, tmp_out);

		lvl10_8[i][0] = tmp_out[0];
		lvl10_8[i][1] = tmp_out[1];
		lvl10_8[i][2] = tmp_out[2];
		lvl10_8[i][3] = tmp_out[3];
		lvl10_8[i][4] = tmp_out[4];
		lvl10_8[i][5] = tmp_out[5];
		lvl10_8[i][6] = tmp_out[6];
		lvl10_8[i][7] = tmp_out[7];

	}

	//11th, 12th, 13th, and 14th levels (8 instantiations of oem16_8)
	for (int i = 0; i < 8; i++)
	{
		t_so tmp_in[16];
		t_so tmp_out[8];

#pragma HLS ARRAY_RESHAPE variable=tmp_in complete dim=0
#pragma HLS ARRAY_RESHAPE variable=tmp_out complete dim=0

		tmp_in[0] = lvl10_8[i * 2][0];
		tmp_in[1] = lvl10_8[i * 2][1];
		tmp_in[2] = lvl10_8[i * 2][2];
		tmp_in[3] = lvl10_8[i * 2][3];
		tmp_in[4] = lvl10_8[i * 2][4];
		tmp_in[5] = lvl10_8[i * 2][5];
		tmp_in[6] = lvl10_8[i * 2][6];
		tmp_in[7] = lvl10_8[i * 2][7];

		tmp_in[8] = lvl10_8[i * 2 + 1][0];
		tmp_in[9] = lvl10_8[i * 2 + 1][1];
		tmp_in[10] = lvl10_8[i * 2 + 1][2];
		tmp_in[11] = lvl10_8[i * 2 + 1][3];
		tmp_in[12] = lvl10_8[i * 2 + 1][4];
		tmp_in[13] = lvl10_8[i * 2 + 1][5];
		tmp_in[14] = lvl10_8[i * 2 + 1][6];
		tmp_in[15] = lvl10_8[i * 2 + 1][7];

		oem16_8(tmp_in, tmp_out);

		lvl14_8[i][0] = tmp_out[0];
		lvl14_8[i][1] = tmp_out[1];
		lvl14_8[i][2] = tmp_out[2];
		lvl14_8[i][3] = tmp_out[3];
		lvl14_8[i][4] = tmp_out[4];
		lvl14_8[i][5] = tmp_out[5];
		lvl14_8[i][6] = tmp_out[6];
		lvl14_8[i][7] = tmp_out[7];
	}

	//15th, 16th, 17th, and 18th levels (4 instantiations of oem16_8)
	for (int i = 0; i < 4; i++)
	{
		t_so tmp_in[16];
		t_so tmp_out[8];

#pragma HLS ARRAY_RESHAPE variable=tmp_in complete dim=0
#pragma HLS ARRAY_RESHAPE variable=tmp_out complete dim=0

		tmp_in[0] = lvl14_8[i * 2][0];
		tmp_in[1] = lvl14_8[i * 2][1];
		tmp_in[2] = lvl14_8[i * 2][2];
		tmp_in[3] = lvl14_8[i * 2][3];
		tmp_in[4] = lvl14_8[i * 2][4];
		tmp_in[5] = lvl14_8[i * 2][5];
		tmp_in[6] = lvl14_8[i * 2][6];
		tmp_in[7] = lvl14_8[i * 2][7];

		tmp_in[8] = lvl14_8[i * 2 + 1][0];
		tmp_in[9] = lvl14_8[i * 2 + 1][1];
		tmp_in[10] = lvl14_8[i * 2 + 1][2];
		tmp_in[11] = lvl14_8[i * 2 + 1][3];
		tmp_in[12] = lvl14_8[i * 2 + 1][4];
		tmp_in[13] = lvl14_8[i * 2 + 1][5];
		tmp_in[14] = lvl14_8[i * 2 + 1][6];
		tmp_in[15] = lvl14_8[i * 2 + 1][7];

		oem16_8(tmp_in, tmp_out);

		lvl18_8[i][0] = tmp_out[0];
		lvl18_8[i][1] = tmp_out[1];
		lvl18_8[i][2] = tmp_out[2];
		lvl18_8[i][3] = tmp_out[3];
		lvl18_8[i][4] = tmp_out[4];
		lvl18_8[i][5] = tmp_out[5];
		lvl18_8[i][6] = tmp_out[6];
		lvl18_8[i][7] = tmp_out[7];
	}

	//19th, 20th, 21st, and 22nd levels (2 instantiations of oem16_8)
	for (int i = 0; i < 2; i++)
	{

		t_so tmp_in[16];
		t_so tmp_out[8];

#pragma HLS ARRAY_RESHAPE variable=tmp_in complete dim=0
#pragma HLS ARRAY_RESHAPE variable=tmp_out complete dim=0

		tmp_in[0] = lvl18_8[i * 2][0];
		tmp_in[1] = lvl18_8[i * 2][1];
		tmp_in[2] = lvl18_8[i * 2][2];
		tmp_in[3] = lvl18_8[i * 2][3];
		tmp_in[4] = lvl18_8[i * 2][4];
		tmp_in[5] = lvl18_8[i * 2][5];
		tmp_in[6] = lvl18_8[i * 2][6];
		tmp_in[7] = lvl18_8[i * 2][7];

		tmp_in[8] = lvl18_8[i * 2 + 1][0];
		tmp_in[9] = lvl18_8[i * 2 + 1][1];
		tmp_in[10] = lvl18_8[i * 2 + 1][2];
		tmp_in[11] = lvl18_8[i * 2 + 1][3];
		tmp_in[12] = lvl18_8[i * 2 + 1][4];
		tmp_in[13] = lvl18_8[i * 2 + 1][5];
		tmp_in[14] = lvl18_8[i * 2 + 1][6];
		tmp_in[15] = lvl18_8[i * 2 + 1][7];

		oem16_8(tmp_in, tmp_out);

		lvl22_8[i][0] = tmp_out[0];
		lvl22_8[i][1] = tmp_out[1];
		lvl22_8[i][2] = tmp_out[2];
		lvl22_8[i][3] = tmp_out[3];
		lvl22_8[i][4] = tmp_out[4];
		lvl22_8[i][5] = tmp_out[5];
		lvl22_8[i][6] = tmp_out[6];
		lvl22_8[i][7] = tmp_out[7];
	}

	lvl22[0] = lvl22_8[0][0];
	lvl22[1] = lvl22_8[0][1];
	lvl22[2] = lvl22_8[0][2];
	lvl22[3] = lvl22_8[0][3];
	lvl22[4] = lvl22_8[0][4];
	lvl22[5] = lvl22_8[0][5];
	lvl22[6] = lvl22_8[0][6];
	lvl22[7] = lvl22_8[0][7];

	lvl22[8] = lvl22_8[1][7];
	lvl22[9] = lvl22_8[1][6];
	lvl22[10] = lvl22_8[1][5];
	lvl22[11] = lvl22_8[1][4];
	lvl22[12] = lvl22_8[1][3];
	lvl22[13] = lvl22_8[1][2];
	lvl22[14] = lvl22_8[1][1];
	lvl22[15] = lvl22_8[1][0];

	max16(lvl22, lvl23);

#if 1
	oem8_f(lvl23, so_out);
#else
	so_out[0] = lvl23[0];
	so_out[1] = lvl23[1];
	so_out[2] = lvl23[2];
	so_out[3] = lvl23[3];
	so_out[4] = lvl23[4];
	so_out[5] = lvl23[5];
	so_out[6] = lvl23[6];
	so_out[7] = lvl23[7];
#endif
	return;
}

void bm16_8_plus(t_so so_in[16], t_so so_out[8])
{
#pragma HLS ARRAY_RESHAPE variable=so_in complete dim=1
#pragma HLS ARRAY_RESHAPE variable=so_out complete dim=1

	t_so lvl0_0[2];
	t_so lvl0_1[2];
	t_so lvl0_2[2];
	t_so lvl0_3[2];
	t_so lvl0_4[2];
	t_so lvl0_5[2];
	t_so lvl0_6[2];
	t_so lvl0_7[2];

	t_so lvl1_0[2];
	t_so lvl1_1[2];
	t_so lvl1_2[2];
	t_so lvl1_3[2];
	t_so lvl1_4[2];
	t_so lvl1_5[2];
	t_so lvl1_6[2];
	t_so lvl1_7[2];

	t_so lvl2[8];

	t_so lvl4[8];

	lvl0_0[0] = so_in[0];
	lvl0_0[1] = so_in[8];

	lvl0_1[0] = so_in[1];
	lvl0_1[1] = so_in[9];

	lvl0_2[0] = so_in[2];
	lvl0_2[1] = so_in[10];

	lvl0_3[0] = so_in[3];
	lvl0_3[1] = so_in[11];

	lvl0_4[0] = so_in[4];
	lvl0_4[1] = so_in[12];

	lvl0_5[0] = so_in[5];
	lvl0_5[1] = so_in[13];

	lvl0_6[0] = so_in[6];
	lvl0_6[1] = so_in[14];

	lvl0_7[0] = so_in[7];
	lvl0_7[1] = so_in[15];

	//first level
	bm2_plus(lvl0_0, lvl1_0);
	bm2_plus(lvl0_1, lvl1_1);
	bm2_plus(lvl0_2, lvl1_2);
	bm2_plus(lvl0_3, lvl1_3);
	bm2_plus(lvl0_4, lvl1_4);
	bm2_plus(lvl0_5, lvl1_5);
	bm2_plus(lvl0_6, lvl1_6);
	bm2_plus(lvl0_7, lvl1_7);

	//second, third and fourth level
	lvl2[0] = lvl1_0[1];
	lvl2[1] = lvl1_1[1];
	lvl2[2] = lvl1_2[1];
	lvl2[3] = lvl1_3[1];
	lvl2[4] = lvl1_4[1];
	lvl2[5] = lvl1_5[1];
	lvl2[6] = lvl1_6[1];
	lvl2[7] = lvl1_7[1];

	bm8_plus(lvl2, lvl4);

	//output assignment
	so_out[0] = lvl4[0];
	so_out[1] = lvl4[1];
	so_out[2] = lvl4[2];
	so_out[3] = lvl4[3];
	so_out[4] = lvl4[4];
	so_out[5] = lvl4[5];
	so_out[6] = lvl4[6];
	so_out[7] = lvl4[7];

	return;
}

void bm2_minus(t_so so_in[2], t_so so_out[2])
{

#pragma HLS ARRAY_RESHAPE variable=so_in complete dim=1
#pragma HLS ARRAY_RESHAPE variable=so_out complete dim=1

	if (so_in[0].et < so_in[1].et)
	{
		so_out[0] = so_in[1];
		so_out[1] = so_in[0];
	}
	else
	{
		so_out[0] = so_in[0];
		so_out[1] = so_in[1];
	}

	return;
}

void bm2_plus(t_so so_in[2], t_so so_out[2])
{
#pragma HLS ARRAY_RESHAPE variable=so_in complete dim=1
#pragma HLS ARRAY_RESHAPE variable=so_out complete dim=1

	if (so_in[0].et < so_in[1].et)
	{
		so_out[0] = so_in[0];
		so_out[1] = so_in[1];
	}
	else
	{
		so_out[0] = so_in[1];
		so_out[1] = so_in[0];
	}

	return;
}

void bm4_minus(t_so so_in[4], t_so so_out[4])
{
#pragma HLS ARRAY_RESHAPE variable=so_in complete dim=1
#pragma HLS ARRAY_RESHAPE variable=so_out complete dim=1

	t_so lvl0_0[2];
	t_so lvl0_1[2];

	t_so lvl1_0[2];
	t_so lvl1_1[2];

	t_so lvl1_0_prim[2];
	t_so lvl1_1_prim[2];

	t_so lvl2_0[2];
	t_so lvl2_1[2];

	lvl0_0[0] = so_in[0];
	lvl0_0[1] = so_in[2];
	lvl0_1[0] = so_in[1];
	lvl0_1[1] = so_in[3];

	//first level
	bm2_minus(lvl0_0, lvl1_0);
	bm2_minus(lvl0_1, lvl1_1);

	lvl1_0_prim[0] = lvl1_0[0];
	lvl1_0_prim[1] = lvl1_1[0];

	lvl1_1_prim[0] = lvl1_0[1];
	lvl1_1_prim[1] = lvl1_1[1];

	//second level
	bm2_minus(lvl1_0_prim, lvl2_0);
	bm2_minus(lvl1_1_prim, lvl2_1);

	//output assignment
	so_out[0] = lvl2_0[0];
	so_out[1] = lvl2_0[1];
	so_out[2] = lvl2_1[0];
	so_out[3] = lvl2_1[1];

	return;
}


void bm4_plus(t_so so_in[4], t_so so_out[4])
{
#pragma HLS ARRAY_RESHAPE variable=so_in complete dim=1
#pragma HLS ARRAY_RESHAPE variable=so_out complete dim=1

	t_so lvl0_0[2];
	t_so lvl0_1[2];

	t_so lvl1_0[2];
	t_so lvl1_1[2];

	t_so lvl1_0_prim[2];
	t_so lvl1_1_prim[2];

	t_so lvl2_0[2];
	t_so lvl2_1[2];

	lvl0_0[0] = so_in[0];
	lvl0_0[1] = so_in[2];
	lvl0_1[0] = so_in[1];
	lvl0_1[1] = so_in[3];

	//first level
	bm2_plus(lvl0_0, lvl1_0);
	bm2_plus(lvl0_1, lvl1_1);

	lvl1_0_prim[0] = lvl1_0[0];
	lvl1_0_prim[1] = lvl1_1[0];

	lvl1_1_prim[0] = lvl1_0[1];
	lvl1_1_prim[1] = lvl1_1[1];

	//second level
	bm2_plus(lvl1_0_prim, lvl2_0);
	bm2_plus(lvl1_1_prim, lvl2_1);

	//output assignment
	so_out[0] = lvl2_0[0];
	so_out[1] = lvl2_0[1];
	so_out[2] = lvl2_1[0];
	so_out[3] = lvl2_1[1];

	return;
}

void bm8_plus(t_so so_in[8], t_so so_out[8])
{
#pragma HLS ARRAY_RESHAPE variable=so_in complete dim=1
#pragma HLS ARRAY_RESHAPE variable=so_out complete dim=1

	t_so lvl0_0[2];
	t_so lvl0_1[2];
	t_so lvl0_2[2];
	t_so lvl0_3[2];

	t_so lvl1_0[2];
	t_so lvl1_1[2];
	t_so lvl1_2[2];
	t_so lvl1_3[2];

	t_so lvl2_0[4];
	t_so lvl2_1[4];

	t_so lvl3_0[4];
	t_so lvl3_1[4];

	lvl0_0[0] = so_in[0];
	lvl0_0[1] = so_in[4];

	lvl0_1[0] = so_in[1];
	lvl0_1[1] = so_in[5];

	lvl0_2[0] = so_in[2];
	lvl0_2[1] = so_in[6];

	lvl0_3[0] = so_in[3];
	lvl0_3[1] = so_in[7];

	//first level
	bm2_plus(lvl0_0, lvl1_0);
	bm2_plus(lvl0_1, lvl1_1);
	bm2_plus(lvl0_2, lvl1_2);
	bm2_plus(lvl0_3, lvl1_3);

	//second and third level
	lvl2_0[0] = lvl1_0[0];
	lvl2_0[1] = lvl1_1[0];
	lvl2_0[2] = lvl1_2[0];
	lvl2_0[3] = lvl1_3[0];

	lvl2_1[0] = lvl1_0[1];
	lvl2_1[1] = lvl1_1[1];
	lvl2_1[2] = lvl1_2[1];
	lvl2_1[3] = lvl1_3[1];

	bm4_plus(lvl2_0, lvl3_0);
	bm4_plus(lvl2_1, lvl3_1);

	//output assignment
	so_out[0] = lvl3_0[0];
	so_out[1] = lvl3_0[1];
	so_out[2] = lvl3_0[2];
	so_out[3] = lvl3_0[3];

	so_out[4] = lvl3_1[0];
	so_out[5] = lvl3_1[1];
	so_out[6] = lvl3_1[2];
	so_out[7] = lvl3_1[3];

	return;
}


void max16(t_so so_in[16], t_so so_out[8])
{
#pragma HLS ARRAY_RESHAPE variable=so_in complete dim=1
#pragma HLS ARRAY_RESHAPE variable=so_out complete dim=1

	t_so lvl0_0[2];
	t_so lvl0_1[2];
	t_so lvl0_2[2];
	t_so lvl0_3[2];
	t_so lvl0_4[2];
	t_so lvl0_5[2];
	t_so lvl0_6[2];
	t_so lvl0_7[2];

	t_so lvl1_0[2];
	t_so lvl1_1[2];
	t_so lvl1_2[2];
	t_so lvl1_3[2];
	t_so lvl1_4[2];
	t_so lvl1_5[2];
	t_so lvl1_6[2];
	t_so lvl1_7[2];

	lvl0_0[0] = so_in[0];
	lvl0_0[1] = so_in[8];

	lvl0_1[0] = so_in[1];
	lvl0_1[1] = so_in[9];

	lvl0_2[0] = so_in[2];
	lvl0_2[1] = so_in[10];

	lvl0_3[0] = so_in[3];
	lvl0_3[1] = so_in[11];

	lvl0_4[0] = so_in[4];
	lvl0_4[1] = so_in[12];

	lvl0_5[0] = so_in[5];
	lvl0_5[1] = so_in[13];

	lvl0_6[0] = so_in[6];
	lvl0_6[1] = so_in[14];

	lvl0_7[0] = so_in[7];
	lvl0_7[1] = so_in[15];

	//first level
	bm2_plus(lvl0_0, lvl1_0);
	bm2_plus(lvl0_1, lvl1_1);
	bm2_plus(lvl0_2, lvl1_2);
	bm2_plus(lvl0_3, lvl1_3);
	bm2_plus(lvl0_4, lvl1_4);
	bm2_plus(lvl0_5, lvl1_5);
	bm2_plus(lvl0_6, lvl1_6);
	bm2_plus(lvl0_7, lvl1_7);


	//output assignment
	so_out[0] = lvl1_0[1];
	so_out[1] = lvl1_1[1];
	so_out[2] = lvl1_2[1];
	so_out[3] = lvl1_3[1];
	so_out[4] = lvl1_4[1];
	so_out[5] = lvl1_5[1];
	so_out[6] = lvl1_6[1];
	so_out[7] = lvl1_7[1];

	return;
}


void oem16_8(t_so so_in[16], t_so so_out[8])
{
#pragma HLS ARRAY_RESHAPE variable=so_in complete dim=1
#pragma HLS ARRAY_RESHAPE variable=so_out complete dim=1

	t_so lvl0[16];

	t_so lvl4[8];

	lvl0[0] = so_in[0];
	lvl0[1] = so_in[1];
	lvl0[2] = so_in[2];
	lvl0[3] = so_in[3];
	lvl0[4] = so_in[4];
	lvl0[5] = so_in[5];
	lvl0[6] = so_in[6];
	lvl0[7] = so_in[7];
	lvl0[8] = so_in[15];
	lvl0[9] = so_in[14];
	lvl0[10] = so_in[13];
	lvl0[11] = so_in[12];
	lvl0[12] = so_in[11];
	lvl0[13] = so_in[10];
	lvl0[14] = so_in[9];
	lvl0[15] = so_in[8];

	//first, second and third level

	bm16_8_plus(lvl0, lvl4);

	//output assignment
	so_out[0] = lvl4[0];
	so_out[1] = lvl4[1];
	so_out[2] = lvl4[2];
	so_out[3] = lvl4[3];
	so_out[4] = lvl4[4];
	so_out[5] = lvl4[5];
	so_out[6] = lvl4[6];
	so_out[7] = lvl4[7];

	return;
}

void oem4(t_so so_in[4], t_so so_out[4])
{
#pragma HLS ARRAY_RESHAPE variable=so_in complete dim=1
#pragma HLS ARRAY_RESHAPE variable=so_out complete dim=1

	t_so lvl0_0[2];
	t_so lvl0_1[2];

	t_so lvl1_0[2];
	t_so lvl1_1[2];

	t_so lvl2[4];

	lvl0_0[0] = so_in[0];
	lvl0_0[1] = so_in[2];
	lvl0_1[0] = so_in[1];
	lvl0_1[1] = so_in[3];

	//first level
	bm2_plus(lvl0_0, lvl1_0);
	bm2_plus(lvl0_1, lvl1_1);

	//second level
	lvl2[0] = lvl1_0[0];
	lvl2[3] = lvl1_1[1];

	lvl2[1] = (lvl1_0[1].et < lvl1_1[0].et ? lvl1_0[1] : lvl1_1[0] );
	lvl2[2] = (lvl1_0[1].et < lvl1_1[0].et ? lvl1_1[0] : lvl1_0[1] );

	//output assignment
	so_out[0] = lvl2[0];
	so_out[1] = lvl2[1];
	so_out[2] = lvl2[2];
	so_out[3] = lvl2[3];

	return;
}

void oem8_f(t_so so_in[8], t_so so_out[8])
{

#pragma HLS ARRAY_RESHAPE variable=so_in complete dim=1
#pragma HLS ARRAY_RESHAPE variable=so_out complete dim=1

	// Level 1
	t_so lvl1_0_in[2];
	t_so lvl1_1_in[2];
	t_so lvl1_2_in[2];
	t_so lvl1_3_in[2];

	t_so lvl1_0_out[2];
	t_so lvl1_1_out[2];
	t_so lvl1_2_out[2];
	t_so lvl1_3_out[2];

	//Level 2
	t_so lvl2_0_in[2];
	t_so lvl2_1_in[2];
	t_so lvl2_2_in[2];
	t_so lvl2_3_in[2];

	t_so lvl2_0_out[2];
	t_so lvl2_1_out[2];
	t_so lvl2_2_out[2];
	t_so lvl2_3_out[2];

	//Level 3
	t_so lvl3_0_in[2];
	t_so lvl3_1_in[2];

	t_so lvl3_0_out[2];
	t_so lvl3_1_out[2];

	//Level 4
	t_so lvl4_0_in[2];
	t_so lvl4_1_in[2];
	t_so lvl4_2_in[2];
	t_so lvl4_3_in[2];

	t_so lvl4_0_out[2];
	t_so lvl4_1_out[2];
	t_so lvl4_2_out[2];
	t_so lvl4_3_out[2];

	//Level 5
	t_so lvl5_0_in[2];
	t_so lvl5_1_in[2];

	t_so lvl5_0_out[2];
	t_so lvl5_1_out[2];

	//Level 6
	t_so lvl6_0_in[2];
	t_so lvl6_1_in[2];
	t_so lvl6_2_in[2];

	t_so lvl6_0_out[2];
	t_so lvl6_1_out[2];
	t_so lvl6_2_out[2];

	//// Level 1
	lvl1_0_in[0]=so_in[0];
	lvl1_0_in[1]=so_in[1];
	lvl1_1_in[0]=so_in[2];
	lvl1_1_in[1]=so_in[3];
	lvl1_2_in[0]=so_in[4];
	lvl1_2_in[1]=so_in[5];
	lvl1_3_in[0]=so_in[6];
	lvl1_3_in[1]=so_in[7];

	bm2_plus(lvl1_0_in, lvl1_0_out);
	bm2_plus(lvl1_1_in, lvl1_1_out);
	bm2_plus(lvl1_2_in, lvl1_2_out);
	bm2_plus(lvl1_3_in, lvl1_3_out);

	//// Level_2
	lvl2_0_in[0]=lvl1_0_out[0];
	lvl2_0_in[1]=lvl1_1_out[0];
	lvl2_1_in[0]=lvl1_0_out[1];
	lvl2_1_in[1]=lvl1_1_out[1];
	lvl2_2_in[0]=lvl1_2_out[0];
	lvl2_2_in[1]=lvl1_3_out[0];
	lvl2_3_in[0]=lvl1_2_out[1];
	lvl2_3_in[1]=lvl1_3_out[1];

	bm2_plus(lvl2_0_in, lvl2_0_out);
	bm2_plus(lvl2_1_in, lvl2_1_out);
	bm2_plus(lvl2_2_in, lvl2_2_out);
	bm2_plus(lvl2_3_in, lvl2_3_out);

	//// Level_3
	lvl3_0_in[0]=lvl2_1_out[0];
	lvl3_0_in[1]=lvl2_0_out[1];
	lvl3_1_in[0]=lvl2_3_out[0];
	lvl3_1_in[1]=lvl2_2_out[1];

	bm2_plus(lvl3_0_in, lvl3_0_out);
	bm2_plus(lvl3_1_in, lvl3_1_out);

	//// Level_4
	lvl4_0_in[0]=lvl2_0_out[0];
	lvl4_0_in[1]=lvl2_2_out[0];
	lvl4_1_in[0]=lvl3_0_out[0];
	lvl4_1_in[1]=lvl3_1_out[0];
	lvl4_2_in[0]=lvl3_0_out[1];
	lvl4_2_in[1]=lvl3_1_out[1];
	lvl4_3_in[0]=lvl2_1_out[1];
	lvl4_3_in[1]=lvl2_3_out[1];

	bm2_plus(lvl4_0_in, lvl4_0_out);
	bm2_plus(lvl4_1_in, lvl4_1_out);
	bm2_plus(lvl4_2_in, lvl4_2_out);
	bm2_plus(lvl4_3_in, lvl4_3_out);

	//// Level_5
	lvl5_0_in[0]=lvl4_2_out[0];
	lvl5_0_in[1]=lvl4_0_out[1];
	lvl5_1_in[0]=lvl4_3_out[0];
	lvl5_1_in[1]=lvl4_1_out[1];

	bm2_plus(lvl5_0_in, lvl5_0_out);
	bm2_plus(lvl5_1_in, lvl5_1_out);

	//// Level_6
	lvl6_0_in[0]=lvl4_1_out[0];
	lvl6_0_in[1]=lvl5_0_out[0];
	lvl6_1_in[0]=lvl5_1_out[0];
	lvl6_1_in[1]=lvl5_0_out[1];
	lvl6_2_in[0]=lvl5_1_out[1];
	lvl6_2_in[1]=lvl4_2_out[1];

	bm2_plus(lvl6_0_in, lvl6_0_out);
	bm2_plus(lvl6_1_in, lvl6_1_out);
	bm2_plus(lvl6_2_in, lvl6_2_out);

	so_out[7] = lvl4_0_out[0];
	so_out[6] = lvl6_0_out[0];
	so_out[5] = lvl6_0_out[1];
	so_out[4] = lvl6_1_out[0];
	so_out[3] = lvl6_1_out[1];
	so_out[2] = lvl6_2_out[0];
	so_out[1] = lvl6_2_out[1];
	so_out[0] = lvl4_3_out[1];

	return;
}

void oem8(t_so so_in[8], t_so so_out[8])
{

#pragma HLS ARRAY_RESHAPE variable=so_in complete dim=1
#pragma HLS ARRAY_RESHAPE variable=so_out complete dim=1

	t_so lvl0_0[4];
	t_so lvl0_1[4];

	t_so lvl2_0[4];
	t_so lvl2_1[4];

	t_so lvl3[8];

	lvl0_0[0] = so_in[0];
	lvl0_0[1] = so_in[2];
	lvl0_0[2] = so_in[4];
	lvl0_0[3] = so_in[6];

	lvl0_1[0] = so_in[1];
	lvl0_1[1] = so_in[3];
	lvl0_1[2] = so_in[5];
	lvl0_1[3] = so_in[7];

	//first and second levels
	oem4(lvl0_0, lvl2_0);
	oem4(lvl0_1, lvl2_1);

	//third level
	lvl3[0] = lvl2_0[0];
	lvl3[7] = lvl2_1[3];

	lvl3[1] = (lvl2_0[1].et < lvl2_1[0].et ? lvl2_0[1] : lvl2_1[0]);
	lvl3[2] = (lvl2_0[1].et < lvl2_1[0].et ? lvl2_1[0] : lvl2_0[1]);

	lvl3[3] = (lvl2_0[2].et < lvl2_1[1].et ? lvl2_0[2] : lvl2_1[1]);
	lvl3[4] = (lvl2_0[2].et < lvl2_1[1].et ? lvl2_1[1] : lvl2_0[2]);

	lvl3[5] = (lvl2_0[3].et < lvl2_1[2].et ? lvl2_0[3] : lvl2_1[2]);
	lvl3[6] = (lvl2_0[3].et < lvl2_1[2].et ? lvl2_1[2] : lvl2_0[3]);

	so_out[0] = lvl3[0];
	so_out[1] = lvl3[1];
	so_out[2] = lvl3[2];
	so_out[3] = lvl3[3];
	so_out[4] = lvl3[4];
	so_out[5] = lvl3[5];
	so_out[6] = lvl3[6];
	so_out[7] = lvl3[7];

	return;
}
