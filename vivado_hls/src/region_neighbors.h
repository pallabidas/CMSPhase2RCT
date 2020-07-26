#ifndef __REGION_NEIGHBORS_H__
#define __REGION_NEIGHBORS_H__

#include "UCTSummaryCard.hpp"

// TODO: consolidate "region_cntr_neighbors_t" table (central regions only)
// with "region_neighbors_t" table (central + forward regions)

typedef struct
{
	// Neighbour coordinates (-1 to indicate no neighbour)
	ap_int<10> nb_C;
	ap_int<10> nb_N;
	ap_int<10> nb_S;
	ap_int<10> nb_E;
	ap_int<10> nb_W;
	ap_int<10> nb_NE;
	ap_int<10> nb_NW;
	ap_int<10> nb_SE;
	ap_int<10> nb_SW;

} region_neighbors_t;

typedef struct {

	// Neighbour coordinates (-1 to indicate no neighbour)
	ap_int<9> nb_reg;

	ap_int<9> nb_N;
	ap_int<9> nb_S;
	ap_int<9> nb_E;
	ap_int<9> nb_W;
	ap_int<9> nb_NE;
	ap_int<9> nb_NW;
	ap_int<9> nb_SE;
	ap_int<9> nb_SW;

} region_cntr_neighbors_t;

const region_neighbors_t rgn_nghbr[NR_CALO_REG] =
{
// Phi:  0
		{ 0,  14, 442,  1, -1, 15, -1, 443,  -1 },
		{ 1,  15, 443,  2,  0, 16, 14, 444, 442 },
		{ 2,  16, 444,  3,  1, 17, 15, 445, 443 },
		{ 3,  17, 445,  4,  2, 18, 16, 446, 444 },
		{ 4,  18, 446,  5,  3, 19, 17, 447, 445 },
		{ 5,  19, 447,  6,  4, 20, 18, 448, 446 },
		{ 6,  20, 448,  7,  5, 21, 19, 449, 447 },
		{ 7,  21, 449,  8,  6, 22, 20, 450, 448 },
		{ 8,  22, 450,  9,  7, 23, 21, 451, 449 },
		{ 9,  23, 451, 10,  8, 24, 22, 452, 450 },
		{ 10, 24, 452, 11,  9, 25, 23, 453, 451 },
		{ 11, 25, 453, 12, 10, 26, 24, 454, 452 },
		{ 12, 26, 454, 13, 11, 27, 25, 455, 453 },
		{ 13, 27, 455, 14, 12, 28, 26, 456, 454 },
		{ 14, 28, 456, 15, 13, 29, 27, 457, 455 },
		{ 15, 29, 457, 16, 14, 30, 28, 458, 456 },
		{ 16, 30, 458, 17, 15, 31, 29, 459, 457 },
		{ 17, 31, 459, 18, 16, 32, 30, 460, 458 },
		{ 18, 32, 460, 19, 17, 33, 31, 461, 459 },
		{ 19, 33, 461, 20, 18, 34, 32, 462, 460 },
		{ 20, 34, 462, 21, 19, 35, 33, 463, 461 },
		{ 21, 35, 463, 22, 20, 36, 34, 464, 462 },
		{ 22, 36, 464, 23, 21, 37, 35, 465, 463 },
		{ 23, 37, 465, 24, 22, 38, 36, 466, 464 },
		{ 24, 38, 466, 25, 23, 39, 37, 467, 465 },
		{ 25, 39, 467, -1, 24, -1, 38, -1, 466 },
		// Phi:  1
		{ 26, 40, 12, 27, -1, 41, -1, 13, -1 },
		{ 27, 41, 13, 28, 26, 42, 40, 14, 12 },
		{ 28, 42, 14, 29, 27, 43, 41, 15, 13 },
		{ 29, 43, 15, 30, 28, 44, 42, 16, 14 },
		{ 30, 44, 16, 31, 29, 45, 43, 17, 15 },
		{ 31, 45, 17, 32, 30, 46, 44, 18, 16 },
		{ 32, 46, 18, 33, 31, 47, 45, 19, 17 },
		{ 33, 47, 19, 34, 32, 48, 46, 20, 18 },
		{ 34, 48, 20, 35, 33, 49, 47, 21, 19 },
		{ 35, 49, 21, 36, 34, 50, 48, 22, 20 },
		{ 36, 50, 22, 37, 35, 51, 49, 23, 21 },
		{ 37, 51, 23, 38, 36, 52, 50, 24, 22 },
		{ 38, 52, 24, 39, 37, 53, 51, 25, 23 },
		{ 39, 53, 25, 40, 38, 54, 52, 26, 24 },
		{ 40, 54, 26, 41, 39, 55, 53, 27, 25 },
		{ 41, 55, 27, 42, 40, 56, 54, 28, 26 },
		{ 42, 56, 28, 43, 41, 57, 55, 29, 27 },
		{ 43, 57, 29, 44, 42, 58, 56, 30, 28 },
		{ 44, 58, 30, 45, 43, 59, 57, 31, 29 },
		{ 45, 59, 31, 46, 44, 60, 58, 32, 30 },
		{ 46, 60, 32, 47, 45, 61, 59, 33, 31 },
		{ 47, 61, 33, 48, 46, 62, 60, 34, 32 },
		{ 48, 62, 34, 49, 47, 63, 61, 35, 33 },
		{ 49, 63, 35, 50, 48, 64, 62, 36, 34 },
		{ 50, 64, 36, 51, 49, 65, 63, 37, 35 },
		{ 51, 65, 37, -1, 50, -1, 64, -1, 36 },
		// Phi:  2
		{ 52, 66, 38, 53, -1, 67, -1, 39, -1 },
		{ 53, 67, 39, 54, 52, 68, 66, 40, 38 },
		{ 54, 68, 40, 55, 53, 69, 67, 41, 39 },
		{ 55, 69, 41, 56, 54, 70, 68, 42, 40 },
		{ 56, 70, 42, 57, 55, 71, 69, 43, 41 },
		{ 57, 71, 43, 58, 56, 72, 70, 44, 42 },
		{ 58, 72, 44, 59, 57, 73, 71, 45, 43 },
		{ 59, 73, 45, 60, 58, 74, 72, 46, 44 },
		{ 60, 74, 46, 61, 59, 75, 73, 47, 45 },
		{ 61, 75, 47, 62, 60, 76, 74, 48, 46 },
		{ 62, 76, 48, 63, 61, 77, 75, 49, 47 },
		{ 63, 77, 49, 64, 62, 78, 76, 50, 48 },
		{ 64, 78, 50, 65, 63, 79, 77, 51, 49 },
		{ 65, 79, 51, 66, 64, 80, 78, 52, 50 },
		{ 66, 80, 52, 67, 65, 81, 79, 53, 51 },
		{ 67, 81, 53, 68, 66, 82, 80, 54, 52 },
		{ 68, 82, 54, 69, 67, 83, 81, 55, 53 },
		{ 69, 83, 55, 70, 68, 84, 82, 56, 54 },
		{ 70, 84, 56, 71, 69, 85, 83, 57, 55 },
		{ 71, 85, 57, 72, 70, 86, 84, 58, 56 },
		{ 72, 86, 58, 73, 71, 87, 85, 59, 57 },
		{ 73, 87, 59, 74, 72, 88, 86, 60, 58 },
		{ 74, 88, 60, 75, 73, 89, 87, 61, 59 },
		{ 75, 89, 61, 76, 74, 90, 88, 62, 60 },
		{ 76, 90, 62, 77, 75, 91, 89, 63, 61 },
		{ 77, 91, 63, -1, 76, -1, 90, -1, 62 },
		// Phi:  3
		{ 78, 92, 64, 79, -1, 93, -1, 65, -1 },
		{ 79, 93, 65, 80, 78, 94, 92, 66, 64 },
		{ 80, 94, 66, 81, 79, 95, 93, 67, 65 },
		{ 81, 95, 67, 82, 80, 96, 94, 68, 66 },
		{ 82, 96, 68, 83, 81, 97, 95, 69, 67 },
		{ 83, 97, 69, 84, 82, 98, 96, 70, 68 },
		{ 84, 98, 70, 85, 83, 99, 97, 71, 69 },
		{ 85, 99, 71, 86, 84, 100, 98, 72, 70 },
		{ 86, 100, 72, 87, 85, 101, 99, 73, 71 },
		{ 87, 101, 73, 88, 86, 102, 100, 74, 72 },
		{ 88, 102, 74, 89, 87, 103, 101, 75, 73 },
		{ 89, 103, 75, 90, 88, 104, 102, 76, 74 },
		{ 90, 104, 76, 91, 89, 105, 103, 77, 75 },
		{ 91, 105, 77, 92, 90, 106, 104, 78, 76 },
		{ 92, 106, 78, 93, 91, 107, 105, 79, 77 },
		{ 93, 107, 79, 94, 92, 108, 106, 80, 78 },
		{ 94, 108, 80, 95, 93, 109, 107, 81, 79 },
		{ 95, 109, 81, 96, 94, 110, 108, 82, 80 },
		{ 96, 110, 82, 97, 95, 111, 109, 83, 81 },
		{ 97, 111, 83, 98, 96, 112, 110, 84, 82 },
		{ 98, 112, 84, 99, 97, 113, 111, 85, 83 },
		{ 99, 113, 85, 100, 98, 114, 112, 86, 84 },
		{ 100, 114, 86, 101, 99, 115, 113, 87, 85 },
		{ 101, 115, 87, 102, 100, 116, 114, 88, 86 },
		{ 102, 116, 88, 103, 101, 117, 115, 89, 87 },
		{ 103, 117, 89, -1, 102, -1, 116, -1, 88 },
		// Phi:  4
		{ 104, 118, 90, 105, -1, 119, -1, 91, -1 },
		{ 105, 119, 91, 106, 104, 120, 118, 92, 90 },
		{ 106, 120, 92, 107, 105, 121, 119, 93, 91 },
		{ 107, 121, 93, 108, 106, 122, 120, 94, 92 },
		{ 108, 122, 94, 109, 107, 123, 121, 95, 93 },
		{ 109, 123, 95, 110, 108, 124, 122, 96, 94 },
		{ 110, 124, 96, 111, 109, 125, 123, 97, 95 },
		{ 111, 125, 97, 112, 110, 126, 124, 98, 96 },
		{ 112, 126, 98, 113, 111, 127, 125, 99, 97 },
		{ 113, 127, 99, 114, 112, 128, 126, 100, 98 },
		{ 114, 128, 100, 115, 113, 129, 127, 101, 99 },
		{ 115, 129, 101, 116, 114, 130, 128, 102, 100 },
		{ 116, 130, 102, 117, 115, 131, 129, 103, 101 },
		{ 117, 131, 103, 118, 116, 132, 130, 104, 102 },
		{ 118, 132, 104, 119, 117, 133, 131, 105, 103 },
		{ 119, 133, 105, 120, 118, 134, 132, 106, 104 },
		{ 120, 134, 106, 121, 119, 135, 133, 107, 105 },
		{ 121, 135, 107, 122, 120, 136, 134, 108, 106 },
		{ 122, 136, 108, 123, 121, 137, 135, 109, 107 },
		{ 123, 137, 109, 124, 122, 138, 136, 110, 108 },
		{ 124, 138, 110, 125, 123, 139, 137, 111, 109 },
		{ 125, 139, 111, 126, 124, 140, 138, 112, 110 },
		{ 126, 140, 112, 127, 125, 141, 139, 113, 111 },
		{ 127, 141, 113, 128, 126, 142, 140, 114, 112 },
		{ 128, 142, 114, 129, 127, 143, 141, 115, 113 },
		{ 129, 143, 115, -1, 128, -1, 142, -1, 114 },
		// Phi:  5
		{ 130, 144, 116, 131, -1, 145, -1, 117, -1 },
		{ 131, 145, 117, 132, 130, 146, 144, 118, 116 },
		{ 132, 146, 118, 133, 131, 147, 145, 119, 117 },
		{ 133, 147, 119, 134, 132, 148, 146, 120, 118 },
		{ 134, 148, 120, 135, 133, 149, 147, 121, 119 },
		{ 135, 149, 121, 136, 134, 150, 148, 122, 120 },
		{ 136, 150, 122, 137, 135, 151, 149, 123, 121 },
		{ 137, 151, 123, 138, 136, 152, 150, 124, 122 },
		{ 138, 152, 124, 139, 137, 153, 151, 125, 123 },
		{ 139, 153, 125, 140, 138, 154, 152, 126, 124 },
		{ 140, 154, 126, 141, 139, 155, 153, 127, 125 },
		{ 141, 155, 127, 142, 140, 156, 154, 128, 126 },
		{ 142, 156, 128, 143, 141, 157, 155, 129, 127 },
		{ 143, 157, 129, 144, 142, 158, 156, 130, 128 },
		{ 144, 158, 130, 145, 143, 159, 157, 131, 129 },
		{ 145, 159, 131, 146, 144, 160, 158, 132, 130 },
		{ 146, 160, 132, 147, 145, 161, 159, 133, 131 },
		{ 147, 161, 133, 148, 146, 162, 160, 134, 132 },
		{ 148, 162, 134, 149, 147, 163, 161, 135, 133 },
		{ 149, 163, 135, 150, 148, 164, 162, 136, 134 },
		{ 150, 164, 136, 151, 149, 165, 163, 137, 135 },
		{ 151, 165, 137, 152, 150, 166, 164, 138, 136 },
		{ 152, 166, 138, 153, 151, 167, 165, 139, 137 },
		{ 153, 167, 139, 154, 152, 168, 166, 140, 138 },
		{ 154, 168, 140, 155, 153, 169, 167, 141, 139 },
		{ 155, 169, 141, -1, 154, -1, 168, -1, 140 },
		// Phi:  6
		{ 156, 170, 142, 157, -1, 171, -1, 143, -1 },
		{ 157, 171, 143, 158, 156, 172, 170, 144, 142 },
		{ 158, 172, 144, 159, 157, 173, 171, 145, 143 },
		{ 159, 173, 145, 160, 158, 174, 172, 146, 144 },
		{ 160, 174, 146, 161, 159, 175, 173, 147, 145 },
		{ 161, 175, 147, 162, 160, 176, 174, 148, 146 },
		{ 162, 176, 148, 163, 161, 177, 175, 149, 147 },
		{ 163, 177, 149, 164, 162, 178, 176, 150, 148 },
		{ 164, 178, 150, 165, 163, 179, 177, 151, 149 },
		{ 165, 179, 151, 166, 164, 180, 178, 152, 150 },
		{ 166, 180, 152, 167, 165, 181, 179, 153, 151 },
		{ 167, 181, 153, 168, 166, 182, 180, 154, 152 },
		{ 168, 182, 154, 169, 167, 183, 181, 155, 153 },
		{ 169, 183, 155, 170, 168, 184, 182, 156, 154 },
		{ 170, 184, 156, 171, 169, 185, 183, 157, 155 },
		{ 171, 185, 157, 172, 170, 186, 184, 158, 156 },
		{ 172, 186, 158, 173, 171, 187, 185, 159, 157 },
		{ 173, 187, 159, 174, 172, 188, 186, 160, 158 },
		{ 174, 188, 160, 175, 173, 189, 187, 161, 159 },
		{ 175, 189, 161, 176, 174, 190, 188, 162, 160 },
		{ 176, 190, 162, 177, 175, 191, 189, 163, 161 },
		{ 177, 191, 163, 178, 176, 192, 190, 164, 162 },
		{ 178, 192, 164, 179, 177, 193, 191, 165, 163 },
		{ 179, 193, 165, 180, 178, 194, 192, 166, 164 },
		{ 180, 194, 166, 181, 179, 195, 193, 167, 165 },
		{ 181, 195, 167, -1, 180, -1, 194, -1, 166 },
		// Phi:  7
		{ 182, 196, 168, 183, -1, 197, -1, 169, -1 },
		{ 183, 197, 169, 184, 182, 198, 196, 170, 168 },
		{ 184, 198, 170, 185, 183, 199, 197, 171, 169 },
		{ 185, 199, 171, 186, 184, 200, 198, 172, 170 },
		{ 186, 200, 172, 187, 185, 201, 199, 173, 171 },
		{ 187, 201, 173, 188, 186, 202, 200, 174, 172 },
		{ 188, 202, 174, 189, 187, 203, 201, 175, 173 },
		{ 189, 203, 175, 190, 188, 204, 202, 176, 174 },
		{ 190, 204, 176, 191, 189, 205, 203, 177, 175 },
		{ 191, 205, 177, 192, 190, 206, 204, 178, 176 },
		{ 192, 206, 178, 193, 191, 207, 205, 179, 177 },
		{ 193, 207, 179, 194, 192, 208, 206, 180, 178 },
		{ 194, 208, 180, 195, 193, 209, 207, 181, 179 },
		{ 195, 209, 181, 196, 194, 210, 208, 182, 180 },
		{ 196, 210, 182, 197, 195, 211, 209, 183, 181 },
		{ 197, 211, 183, 198, 196, 212, 210, 184, 182 },
		{ 198, 212, 184, 199, 197, 213, 211, 185, 183 },
		{ 199, 213, 185, 200, 198, 214, 212, 186, 184 },
		{ 200, 214, 186, 201, 199, 215, 213, 187, 185 },
		{ 201, 215, 187, 202, 200, 216, 214, 188, 186 },
		{ 202, 216, 188, 203, 201, 217, 215, 189, 187 },
		{ 203, 217, 189, 204, 202, 218, 216, 190, 188 },
		{ 204, 218, 190, 205, 203, 219, 217, 191, 189 },
		{ 205, 219, 191, 206, 204, 220, 218, 192, 190 },
		{ 206, 220, 192, 207, 205, 221, 219, 193, 191 },
		{ 207, 221, 193, -1, 206, -1, 220, -1, 192 },
		// Phi:  8
		{ 208, 222, 194, 209, -1, 223, -1, 195, -1 },
		{ 209, 223, 195, 210, 208, 224, 222, 196, 194 },
		{ 210, 224, 196, 211, 209, 225, 223, 197, 195 },
		{ 211, 225, 197, 212, 210, 226, 224, 198, 196 },
		{ 212, 226, 198, 213, 211, 227, 225, 199, 197 },
		{ 213, 227, 199, 214, 212, 228, 226, 200, 198 },
		{ 214, 228, 200, 215, 213, 229, 227, 201, 199 },
		{ 215, 229, 201, 216, 214, 230, 228, 202, 200 },
		{ 216, 230, 202, 217, 215, 231, 229, 203, 201 },
		{ 217, 231, 203, 218, 216, 232, 230, 204, 202 },
		{ 218, 232, 204, 219, 217, 233, 231, 205, 203 },
		{ 219, 233, 205, 220, 218, 234, 232, 206, 204 },
		{ 220, 234, 206, 221, 219, 235, 233, 207, 205 },
		{ 221, 235, 207, 222, 220, 236, 234, 208, 206 },
		{ 222, 236, 208, 223, 221, 237, 235, 209, 207 },
		{ 223, 237, 209, 224, 222, 238, 236, 210, 208 },
		{ 224, 238, 210, 225, 223, 239, 237, 211, 209 },
		{ 225, 239, 211, 226, 224, 240, 238, 212, 210 },
		{ 226, 240, 212, 227, 225, 241, 239, 213, 211 },
		{ 227, 241, 213, 228, 226, 242, 240, 214, 212 },
		{ 228, 242, 214, 229, 227, 243, 241, 215, 213 },
		{ 229, 243, 215, 230, 228, 244, 242, 216, 214 },
		{ 230, 244, 216, 231, 229, 245, 243, 217, 215 },
		{ 231, 245, 217, 232, 230, 246, 244, 218, 216 },
		{ 232, 246, 218, 233, 231, 247, 245, 219, 217 },
		{ 233, 247, 219, -1, 232, -1, 246, -1, 218 },
		// Phi:  9
		{ 234, 248, 220, 235, -1, 249, -1, 221, -1 },
		{ 235, 249, 221, 236, 234, 250, 248, 222, 220 },
		{ 236, 250, 222, 237, 235, 251, 249, 223, 221 },
		{ 237, 251, 223, 238, 236, 252, 250, 224, 222 },
		{ 238, 252, 224, 239, 237, 253, 251, 225, 223 },
		{ 239, 253, 225, 240, 238, 254, 252, 226, 224 },
		{ 240, 254, 226, 241, 239, 255, 253, 227, 225 },
		{ 241, 255, 227, 242, 240, 256, 254, 228, 226 },
		{ 242, 256, 228, 243, 241, 257, 255, 229, 227 },
		{ 243, 257, 229, 244, 242, 258, 256, 230, 228 },
		{ 244, 258, 230, 245, 243, 259, 257, 231, 229 },
		{ 245, 259, 231, 246, 244, 260, 258, 232, 230 },
		{ 246, 260, 232, 247, 245, 261, 259, 233, 231 },
		{ 247, 261, 233, 248, 246, 262, 260, 234, 232 },
		{ 248, 262, 234, 249, 247, 263, 261, 235, 233 },
		{ 249, 263, 235, 250, 248, 264, 262, 236, 234 },
		{ 250, 264, 236, 251, 249, 265, 263, 237, 235 },
		{ 251, 265, 237, 252, 250, 266, 264, 238, 236 },
		{ 252, 266, 238, 253, 251, 267, 265, 239, 237 },
		{ 253, 267, 239, 254, 252, 268, 266, 240, 238 },
		{ 254, 268, 240, 255, 253, 269, 267, 241, 239 },
		{ 255, 269, 241, 256, 254, 270, 268, 242, 240 },
		{ 256, 270, 242, 257, 255, 271, 269, 243, 241 },
		{ 257, 271, 243, 258, 256, 272, 270, 244, 242 },
		{ 258, 272, 244, 259, 257, 273, 271, 245, 243 },
		{ 259, 273, 245, -1, 258, -1, 272, -1, 244 },
		// Phi:  10
		{ 260, 274, 246, 261, -1, 275, -1, 247, -1 },
		{ 261, 275, 247, 262, 260, 276, 274, 248, 246 },
		{ 262, 276, 248, 263, 261, 277, 275, 249, 247 },
		{ 263, 277, 249, 264, 262, 278, 276, 250, 248 },
		{ 264, 278, 250, 265, 263, 279, 277, 251, 249 },
		{ 265, 279, 251, 266, 264, 280, 278, 252, 250 },
		{ 266, 280, 252, 267, 265, 281, 279, 253, 251 },
		{ 267, 281, 253, 268, 266, 282, 280, 254, 252 },
		{ 268, 282, 254, 269, 267, 283, 281, 255, 253 },
		{ 269, 283, 255, 270, 268, 284, 282, 256, 254 },
		{ 270, 284, 256, 271, 269, 285, 283, 257, 255 },
		{ 271, 285, 257, 272, 270, 286, 284, 258, 256 },
		{ 272, 286, 258, 273, 271, 287, 285, 259, 257 },
		{ 273, 287, 259, 274, 272, 288, 286, 260, 258 },
		{ 274, 288, 260, 275, 273, 289, 287, 261, 259 },
		{ 275, 289, 261, 276, 274, 290, 288, 262, 260 },
		{ 276, 290, 262, 277, 275, 291, 289, 263, 261 },
		{ 277, 291, 263, 278, 276, 292, 290, 264, 262 },
		{ 278, 292, 264, 279, 277, 293, 291, 265, 263 },
		{ 279, 293, 265, 280, 278, 294, 292, 266, 264 },
		{ 280, 294, 266, 281, 279, 295, 293, 267, 265 },
		{ 281, 295, 267, 282, 280, 296, 294, 268, 266 },
		{ 282, 296, 268, 283, 281, 297, 295, 269, 267 },
		{ 283, 297, 269, 284, 282, 298, 296, 270, 268 },
		{ 284, 298, 270, 285, 283, 299, 297, 271, 269 },
		{ 285, 299, 271, -1, 284, -1, 298, -1, 270 },
		// Phi:  11
		{ 286, 300, 272, 287, -1, 301, -1, 273, -1 },
		{ 287, 301, 273, 288, 286, 302, 300, 274, 272 },
		{ 288, 302, 274, 289, 287, 303, 301, 275, 273 },
		{ 289, 303, 275, 290, 288, 304, 302, 276, 274 },
		{ 290, 304, 276, 291, 289, 305, 303, 277, 275 },
		{ 291, 305, 277, 292, 290, 306, 304, 278, 276 },
		{ 292, 306, 278, 293, 291, 307, 305, 279, 277 },
		{ 293, 307, 279, 294, 292, 308, 306, 280, 278 },
		{ 294, 308, 280, 295, 293, 309, 307, 281, 279 },
		{ 295, 309, 281, 296, 294, 310, 308, 282, 280 },
		{ 296, 310, 282, 297, 295, 311, 309, 283, 281 },
		{ 297, 311, 283, 298, 296, 312, 310, 284, 282 },
		{ 298, 312, 284, 299, 297, 313, 311, 285, 283 },
		{ 299, 313, 285, 300, 298, 314, 312, 286, 284 },
		{ 300, 314, 286, 301, 299, 315, 313, 287, 285 },
		{ 301, 315, 287, 302, 300, 316, 314, 288, 286 },
		{ 302, 316, 288, 303, 301, 317, 315, 289, 287 },
		{ 303, 317, 289, 304, 302, 318, 316, 290, 288 },
		{ 304, 318, 290, 305, 303, 319, 317, 291, 289 },
		{ 305, 319, 291, 306, 304, 320, 318, 292, 290 },
		{ 306, 320, 292, 307, 305, 321, 319, 293, 291 },
		{ 307, 321, 293, 308, 306, 322, 320, 294, 292 },
		{ 308, 322, 294, 309, 307, 323, 321, 295, 293 },
		{ 309, 323, 295, 310, 308, 324, 322, 296, 294 },
		{ 310, 324, 296, 311, 309, 325, 323, 297, 295 },
		{ 311, 325, 297, -1, 310, -1, 324, -1, 296 },
		// Phi:  12
		{ 312, 326, 298, 313, -1, 327, -1, 299, -1 },
		{ 313, 327, 299, 314, 312, 328, 326, 300, 298 },
		{ 314, 328, 300, 315, 313, 329, 327, 301, 299 },
		{ 315, 329, 301, 316, 314, 330, 328, 302, 300 },
		{ 316, 330, 302, 317, 315, 331, 329, 303, 301 },
		{ 317, 331, 303, 318, 316, 332, 330, 304, 302 },
		{ 318, 332, 304, 319, 317, 333, 331, 305, 303 },
		{ 319, 333, 305, 320, 318, 334, 332, 306, 304 },
		{ 320, 334, 306, 321, 319, 335, 333, 307, 305 },
		{ 321, 335, 307, 322, 320, 336, 334, 308, 306 },
		{ 322, 336, 308, 323, 321, 337, 335, 309, 307 },
		{ 323, 337, 309, 324, 322, 338, 336, 310, 308 },
		{ 324, 338, 310, 325, 323, 339, 337, 311, 309 },
		{ 325, 339, 311, 326, 324, 340, 338, 312, 310 },
		{ 326, 340, 312, 327, 325, 341, 339, 313, 311 },
		{ 327, 341, 313, 328, 326, 342, 340, 314, 312 },
		{ 328, 342, 314, 329, 327, 343, 341, 315, 313 },
		{ 329, 343, 315, 330, 328, 344, 342, 316, 314 },
		{ 330, 344, 316, 331, 329, 345, 343, 317, 315 },
		{ 331, 345, 317, 332, 330, 346, 344, 318, 316 },
		{ 332, 346, 318, 333, 331, 347, 345, 319, 317 },
		{ 333, 347, 319, 334, 332, 348, 346, 320, 318 },
		{ 334, 348, 320, 335, 333, 349, 347, 321, 319 },
		{ 335, 349, 321, 336, 334, 350, 348, 322, 320 },
		{ 336, 350, 322, 337, 335, 351, 349, 323, 321 },
		{ 337, 351, 323, -1, 336, -1, 350, -1, 322 },
		// Phi:  13
		{ 338, 352, 324, 339, -1, 353, -1, 325, -1 },
		{ 339, 353, 325, 340, 338, 354, 352, 326, 324 },
		{ 340, 354, 326, 341, 339, 355, 353, 327, 325 },
		{ 341, 355, 327, 342, 340, 356, 354, 328, 326 },
		{ 342, 356, 328, 343, 341, 357, 355, 329, 327 },
		{ 343, 357, 329, 344, 342, 358, 356, 330, 328 },
		{ 344, 358, 330, 345, 343, 359, 357, 331, 329 },
		{ 345, 359, 331, 346, 344, 360, 358, 332, 330 },
		{ 346, 360, 332, 347, 345, 361, 359, 333, 331 },
		{ 347, 361, 333, 348, 346, 362, 360, 334, 332 },
		{ 348, 362, 334, 349, 347, 363, 361, 335, 333 },
		{ 349, 363, 335, 350, 348, 364, 362, 336, 334 },
		{ 350, 364, 336, 351, 349, 365, 363, 337, 335 },
		{ 351, 365, 337, 352, 350, 366, 364, 338, 336 },
		{ 352, 366, 338, 353, 351, 367, 365, 339, 337 },
		{ 353, 367, 339, 354, 352, 368, 366, 340, 338 },
		{ 354, 368, 340, 355, 353, 369, 367, 341, 339 },
		{ 355, 369, 341, 356, 354, 370, 368, 342, 340 },
		{ 356, 370, 342, 357, 355, 371, 369, 343, 341 },
		{ 357, 371, 343, 358, 356, 372, 370, 344, 342 },
		{ 358, 372, 344, 359, 357, 373, 371, 345, 343 },
		{ 359, 373, 345, 360, 358, 374, 372, 346, 344 },
		{ 360, 374, 346, 361, 359, 375, 373, 347, 345 },
		{ 361, 375, 347, 362, 360, 376, 374, 348, 346 },
		{ 362, 376, 348, 363, 361, 377, 375, 349, 347 },
		{ 363, 377, 349, -1, 362, -1, 376, -1, 348 },
		// Phi:  14
		{ 364, 378, 350, 365, -1, 379, -1, 351, -1 },
		{ 365, 379, 351, 366, 364, 380, 378, 352, 350 },
		{ 366, 380, 352, 367, 365, 381, 379, 353, 351 },
		{ 367, 381, 353, 368, 366, 382, 380, 354, 352 },
		{ 368, 382, 354, 369, 367, 383, 381, 355, 353 },
		{ 369, 383, 355, 370, 368, 384, 382, 356, 354 },
		{ 370, 384, 356, 371, 369, 385, 383, 357, 355 },
		{ 371, 385, 357, 372, 370, 386, 384, 358, 356 },
		{ 372, 386, 358, 373, 371, 387, 385, 359, 357 },
		{ 373, 387, 359, 374, 372, 388, 386, 360, 358 },
		{ 374, 388, 360, 375, 373, 389, 387, 361, 359 },
		{ 375, 389, 361, 376, 374, 390, 388, 362, 360 },
		{ 376, 390, 362, 377, 375, 391, 389, 363, 361 },
		{ 377, 391, 363, 378, 376, 392, 390, 364, 362 },
		{ 378, 392, 364, 379, 377, 393, 391, 365, 363 },
		{ 379, 393, 365, 380, 378, 394, 392, 366, 364 },
		{ 380, 394, 366, 381, 379, 395, 393, 367, 365 },
		{ 381, 395, 367, 382, 380, 396, 394, 368, 366 },
		{ 382, 396, 368, 383, 381, 397, 395, 369, 367 },
		{ 383, 397, 369, 384, 382, 398, 396, 370, 368 },
		{ 384, 398, 370, 385, 383, 399, 397, 371, 369 },
		{ 385, 399, 371, 386, 384, 400, 398, 372, 370 },
		{ 386, 400, 372, 387, 385, 401, 399, 373, 371 },
		{ 387, 401, 373, 388, 386, 402, 400, 374, 372 },
		{ 388, 402, 374, 389, 387, 403, 401, 375, 373 },
		{ 389, 403, 375, -1, 388, -1, 402, -1, 374 },
		// Phi:  15
		{ 390, 404, 376, 391, -1, 405, -1, 377, -1 },
		{ 391, 405, 377, 392, 390, 406, 404, 378, 376 },
		{ 392, 406, 378, 393, 391, 407, 405, 379, 377 },
		{ 393, 407, 379, 394, 392, 408, 406, 380, 378 },
		{ 394, 408, 380, 395, 393, 409, 407, 381, 379 },
		{ 395, 409, 381, 396, 394, 410, 408, 382, 380 },
		{ 396, 410, 382, 397, 395, 411, 409, 383, 381 },
		{ 397, 411, 383, 398, 396, 412, 410, 384, 382 },
		{ 398, 412, 384, 399, 397, 413, 411, 385, 383 },
		{ 399, 413, 385, 400, 398, 414, 412, 386, 384 },
		{ 400, 414, 386, 401, 399, 415, 413, 387, 385 },
		{ 401, 415, 387, 402, 400, 416, 414, 388, 386 },
		{ 402, 416, 388, 403, 401, 417, 415, 389, 387 },
		{ 403, 417, 389, 404, 402, 418, 416, 390, 388 },
		{ 404, 418, 390, 405, 403, 419, 417, 391, 389 },
		{ 405, 419, 391, 406, 404, 420, 418, 392, 390 },
		{ 406, 420, 392, 407, 405, 421, 419, 393, 391 },
		{ 407, 421, 393, 408, 406, 422, 420, 394, 392 },
		{ 408, 422, 394, 409, 407, 423, 421, 395, 393 },
		{ 409, 423, 395, 410, 408, 424, 422, 396, 394 },
		{ 410, 424, 396, 411, 409, 425, 423, 397, 395 },
		{ 411, 425, 397, 412, 410, 426, 424, 398, 396 },
		{ 412, 426, 398, 413, 411, 427, 425, 399, 397 },
		{ 413, 427, 399, 414, 412, 428, 426, 400, 398 },
		{ 414, 428, 400, 415, 413, 429, 427, 401, 399 },
		{ 415, 429, 401, -1, 414, -1, 428, -1, 400 },
		// Phi:  16
		{ 416, 430, 402, 417, -1, 431, -1, 403, -1 },
		{ 417, 431, 403, 418, 416, 432, 430, 404, 402 },
		{ 418, 432, 404, 419, 417, 433, 431, 405, 403 },
		{ 419, 433, 405, 420, 418, 434, 432, 406, 404 },
		{ 420, 434, 406, 421, 419, 435, 433, 407, 405 },
		{ 421, 435, 407, 422, 420, 436, 434, 408, 406 },
		{ 422, 436, 408, 423, 421, 437, 435, 409, 407 },
		{ 423, 437, 409, 424, 422, 438, 436, 410, 408 },
		{ 424, 438, 410, 425, 423, 439, 437, 411, 409 },
		{ 425, 439, 411, 426, 424, 440, 438, 412, 410 },
		{ 426, 440, 412, 427, 425, 441, 439, 413, 411 },
		{ 427, 441, 413, 428, 426, 442, 440, 414, 412 },
		{ 428, 442, 414, 429, 427, 443, 441, 415, 413 },
		{ 429, 443, 415, 430, 428, 444, 442, 416, 414 },
		{ 430, 444, 416, 431, 429, 445, 443, 417, 415 },
		{ 431, 445, 417, 432, 430, 446, 444, 418, 416 },
		{ 432, 446, 418, 433, 431, 447, 445, 419, 417 },
		{ 433, 447, 419, 434, 432, 448, 446, 420, 418 },
		{ 434, 448, 420, 435, 433, 449, 447, 421, 419 },
		{ 435, 449, 421, 436, 434, 450, 448, 422, 420 },
		{ 436, 450, 422, 437, 435, 451, 449, 423, 421 },
		{ 437, 451, 423, 438, 436, 452, 450, 424, 422 },
		{ 438, 452, 424, 439, 437, 453, 451, 425, 423 },
		{ 439, 453, 425, 440, 438, 454, 452, 426, 424 },
		{ 440, 454, 426, 441, 439, 455, 453, 427, 425 },
		{ 441, 455, 427, -1, 440, -1, 454, -1, 426 },
		// Phi:  17
		{ 442, 0, 428, 443, -1, 1, -1, 429, -1 },
		{ 443, 1, 429, 444, 442, 2, 0, 430, 428 },
		{ 444, 2, 430, 445, 443, 3, 1, 431, 429 },
		{ 445, 3, 431, 446, 444, 4, 2, 432, 430 },
		{ 446, 4, 432, 447, 445, 5, 3, 433, 431 },
		{ 447, 5, 433, 448, 446, 6, 4, 434, 432 },
		{ 448, 6, 434, 449, 447, 7, 5, 435, 433 },
		{ 449, 7, 435, 450, 448, 8, 6, 436, 434 },
		{ 450, 8, 436, 451, 449, 9, 7, 437, 435 },
		{ 451, 9, 437, 452, 450, 10, 8, 438, 436 },
		{ 452, 10, 438, 453, 451, 11, 9, 439, 437 },
		{ 453, 11, 439, 454, 452, 12, 10, 440, 438 },
		{ 454, 12, 440, 455, 453, 13, 11, 441, 439 },
		{ 455, 13, 441, 456, 454, 14, 12, 442, 440 },
		{ 456, 14, 442, 457, 455, 15, 13, 443, 441 },
		{ 457, 15, 443, 458, 456, 16, 14, 444, 442 },
		{ 458, 16, 444, 459, 457, 17, 15, 445, 443 },
		{ 459, 17, 445, 460, 458, 18, 16, 446, 444 },
		{ 460, 18, 446, 461, 459, 19, 17, 447, 445 },
		{ 461, 19, 447, 462, 460, 20, 18, 448, 446 },
		{ 462, 20, 448, 463, 461, 21, 19, 449, 447 },
		{ 463, 21, 449, 464, 462, 22, 20, 450, 448 },
		{ 464, 22, 450, 465, 463, 23, 21, 451, 449 },
		{ 465, 23, 451, 466, 464, 24, 22, 452, 450 },
		{ 466, 24, 452, 467, 465, 25, 23, 453, 451 },
		{ 467, 25, 453, -1, 466, -1, 24, -1, 452 } };

////



const region_cntr_neighbors_t region_cntr_neighbors[252] = {

		{ 0 , 14 , 238 , 1 , -1 , 15 , -1 , 239 , -1 },
		{ 1 , 15 , 239 , 2 , 0 , 16 , 14 , 240 , 238 },
		{ 2 , 16 , 240 , 3 , 1 , 17 , 15 , 241 , 239 },
		{ 3 , 17 , 241 , 4 , 2 , 18 , 16 , 242 , 240 },
		{ 4 , 18 , 242 , 5 , 3 , 19 , 17 , 243 , 241 },
		{ 5 , 19 , 243 , 6 , 4 , 20 , 18 , 244 , 242 },
		{ 6 , 20 , 244 , 7 , 5 , 21 , 19 , 245 , 243 },
		{ 7 , 21 , 245 , 8 , 6 , 22 , 20 , 246 , 244 },
		{ 8 , 22 , 246 , 9 , 7 , 23 , 21 , 247 , 245 },
		{ 9 , 23 , 247 , 10 , 8 , 24 , 22 , 248 , 246 },
		{ 10 , 24 , 248 , 11 , 9 , 25 , 23 , 249 , 247 },
		{ 11 , 25 , 249 , 12 , 10 , 26 , 24 , 250 , 248 },
		{ 12 , 26 , 250 , 13 , 11 , 27 , 25 , 251 , 249 },
		{ 13 , 27 , 251 , -1 , 12 , -1 , 26 , -1 , 250 },
		{ 14 , 28 , 0 , 15 , -1 , 29 , -1 , 1 , -1 },
		{ 15 , 29 , 1 , 16 , 14 , 30 , 28 , 2 , 0 },
		{ 16 , 30 , 2 , 17 , 15 , 31 , 29 , 3 , 1 },
		{ 17 , 31 , 3 , 18 , 16 , 32 , 30 , 4 , 2 },
		{ 18 , 32 , 4 , 19 , 17 , 33 , 31 , 5 , 3 },
		{ 19 , 33 , 5 , 20 , 18 , 34 , 32 , 6 , 4 },
		{ 20 , 34 , 6 , 21 , 19 , 35 , 33 , 7 , 5 },
		{ 21 , 35 , 7 , 22 , 20 , 36 , 34 , 8 , 6 },
		{ 22 , 36 , 8 , 23 , 21 , 37 , 35 , 9 , 7 },
		{ 23 , 37 , 9 , 24 , 22 , 38 , 36 , 10 , 8 },
		{ 24 , 38 , 10 , 25 , 23 , 39 , 37 , 11 , 9 },
		{ 25 , 39 , 11 , 26 , 24 , 40 , 38 , 12 , 10 },
		{ 26 , 40 , 12 , 27 , 25 , 41 , 39 , 13 , 11 },
		{ 27 , 41 , 13 , -1 , 26 , -1 , 40 , -1 , 12 },
		{ 28 , 42 , 14 , 29 , -1 , 43 , -1 , 15 , -1 },
		{ 29 , 43 , 15 , 30 , 28 , 44 , 42 , 16 , 14 },
		{ 30 , 44 , 16 , 31 , 29 , 45 , 43 , 17 , 15 },
		{ 31 , 45 , 17 , 32 , 30 , 46 , 44 , 18 , 16 },
		{ 32 , 46 , 18 , 33 , 31 , 47 , 45 , 19 , 17 },
		{ 33 , 47 , 19 , 34 , 32 , 48 , 46 , 20 , 18 },
		{ 34 , 48 , 20 , 35 , 33 , 49 , 47 , 21 , 19 },
		{ 35 , 49 , 21 , 36 , 34 , 50 , 48 , 22 , 20 },
		{ 36 , 50 , 22 , 37 , 35 , 51 , 49 , 23 , 21 },
		{ 37 , 51 , 23 , 38 , 36 , 52 , 50 , 24 , 22 },
		{ 38 , 52 , 24 , 39 , 37 , 53 , 51 , 25 , 23 },
		{ 39 , 53 , 25 , 40 , 38 , 54 , 52 , 26 , 24 },
		{ 40 , 54 , 26 , 41 , 39 , 55 , 53 , 27 , 25 },
		{ 41 , 55 , 27 , -1 , 40 , -1 , 54 , -1 , 26 },
		{ 42 , 56 , 28 , 43 , -1 , 57 , -1 , 29 , -1 },
		{ 43 , 57 , 29 , 44 , 42 , 58 , 56 , 30 , 28 },
		{ 44 , 58 , 30 , 45 , 43 , 59 , 57 , 31 , 29 },
		{ 45 , 59 , 31 , 46 , 44 , 60 , 58 , 32 , 30 },
		{ 46 , 60 , 32 , 47 , 45 , 61 , 59 , 33 , 31 },
		{ 47 , 61 , 33 , 48 , 46 , 62 , 60 , 34 , 32 },
		{ 48 , 62 , 34 , 49 , 47 , 63 , 61 , 35 , 33 },
		{ 49 , 63 , 35 , 50 , 48 , 64 , 62 , 36 , 34 },
		{ 50 , 64 , 36 , 51 , 49 , 65 , 63 , 37 , 35 },
		{ 51 , 65 , 37 , 52 , 50 , 66 , 64 , 38 , 36 },
		{ 52 , 66 , 38 , 53 , 51 , 67 , 65 , 39 , 37 },
		{ 53 , 67 , 39 , 54 , 52 , 68 , 66 , 40 , 38 },
		{ 54 , 68 , 40 , 55 , 53 , 69 , 67 , 41 , 39 },
		{ 55 , 69 , 41 , -1 , 54 , -1 , 68 , -1 , 40 },
		{ 56 , 70 , 42 , 57 , -1 , 71 , -1 , 43 , -1 },
		{ 57 , 71 , 43 , 58 , 56 , 72 , 70 , 44 , 42 },
		{ 58 , 72 , 44 , 59 , 57 , 73 , 71 , 45 , 43 },
		{ 59 , 73 , 45 , 60 , 58 , 74 , 72 , 46 , 44 },
		{ 60 , 74 , 46 , 61 , 59 , 75 , 73 , 47 , 45 },
		{ 61 , 75 , 47 , 62 , 60 , 76 , 74 , 48 , 46 },
		{ 62 , 76 , 48 , 63 , 61 , 77 , 75 , 49 , 47 },
		{ 63 , 77 , 49 , 64 , 62 , 78 , 76 , 50 , 48 },
		{ 64 , 78 , 50 , 65 , 63 , 79 , 77 , 51 , 49 },
		{ 65 , 79 , 51 , 66 , 64 , 80 , 78 , 52 , 50 },
		{ 66 , 80 , 52 , 67 , 65 , 81 , 79 , 53 , 51 },
		{ 67 , 81 , 53 , 68 , 66 , 82 , 80 , 54 , 52 },
		{ 68 , 82 , 54 , 69 , 67 , 83 , 81 , 55 , 53 },
		{ 69 , 83 , 55 , -1 , 68 , -1 , 82 , -1 , 54 },
		{ 70 , 84 , 56 , 71 , -1 , 85 , -1 , 57 , -1 },
		{ 71 , 85 , 57 , 72 , 70 , 86 , 84 , 58 , 56 },
		{ 72 , 86 , 58 , 73 , 71 , 87 , 85 , 59 , 57 },
		{ 73 , 87 , 59 , 74 , 72 , 88 , 86 , 60 , 58 },
		{ 74 , 88 , 60 , 75 , 73 , 89 , 87 , 61 , 59 },
		{ 75 , 89 , 61 , 76 , 74 , 90 , 88 , 62 , 60 },
		{ 76 , 90 , 62 , 77 , 75 , 91 , 89 , 63 , 61 },
		{ 77 , 91 , 63 , 78 , 76 , 92 , 90 , 64 , 62 },
		{ 78 , 92 , 64 , 79 , 77 , 93 , 91 , 65 , 63 },
		{ 79 , 93 , 65 , 80 , 78 , 94 , 92 , 66 , 64 },
		{ 80 , 94 , 66 , 81 , 79 , 95 , 93 , 67 , 65 },
		{ 81 , 95 , 67 , 82 , 80 , 96 , 94 , 68 , 66 },
		{ 82 , 96 , 68 , 83 , 81 , 97 , 95 , 69 , 67 },
		{ 83 , 97 , 69 , -1 , 82 , -1 , 96 , -1 , 68 },
		{ 84 , 98 , 70 , 85 , -1 , 99 , -1 , 71 , -1 },
		{ 85 , 99 , 71 , 86 , 84 , 100 , 98 , 72 , 70 },
		{ 86 , 100 , 72 , 87 , 85 , 101 , 99 , 73 , 71 },
		{ 87 , 101 , 73 , 88 , 86 , 102 , 100 , 74 , 72 },
		{ 88 , 102 , 74 , 89 , 87 , 103 , 101 , 75 , 73 },
		{ 89 , 103 , 75 , 90 , 88 , 104 , 102 , 76 , 74 },
		{ 90 , 104 , 76 , 91 , 89 , 105 , 103 , 77 , 75 },
		{ 91 , 105 , 77 , 92 , 90 , 106 , 104 , 78 , 76 },
		{ 92 , 106 , 78 , 93 , 91 , 107 , 105 , 79 , 77 },
		{ 93 , 107 , 79 , 94 , 92 , 108 , 106 , 80 , 78 },
		{ 94 , 108 , 80 , 95 , 93 , 109 , 107 , 81 , 79 },
		{ 95 , 109 , 81 , 96 , 94 , 110 , 108 , 82 , 80 },
		{ 96 , 110 , 82 , 97 , 95 , 111 , 109 , 83 , 81 },
		{ 97 , 111 , 83 , -1 , 96 , -1 , 110 , -1 , 82 },
		{ 98 , 112 , 84 , 99 , -1 , 113 , -1 , 85 , -1 },
		{ 99 , 113 , 85 , 100 , 98 , 114 , 112 , 86 , 84 },
		{ 100 , 114 , 86 , 101 , 99 , 115 , 113 , 87 , 85 },
		{ 101 , 115 , 87 , 102 , 100 , 116 , 114 , 88 , 86 },
		{ 102 , 116 , 88 , 103 , 101 , 117 , 115 , 89 , 87 },
		{ 103 , 117 , 89 , 104 , 102 , 118 , 116 , 90 , 88 },
		{ 104 , 118 , 90 , 105 , 103 , 119 , 117 , 91 , 89 },
		{ 105 , 119 , 91 , 106 , 104 , 120 , 118 , 92 , 90 },
		{ 106 , 120 , 92 , 107 , 105 , 121 , 119 , 93 , 91 },
		{ 107 , 121 , 93 , 108 , 106 , 122 , 120 , 94 , 92 },
		{ 108 , 122 , 94 , 109 , 107 , 123 , 121 , 95 , 93 },
		{ 109 , 123 , 95 , 110 , 108 , 124 , 122 , 96 , 94 },
		{ 110 , 124 , 96 , 111 , 109 , 125 , 123 , 97 , 95 },
		{ 111 , 125 , 97 , -1 , 110 , -1 , 124 , -1 , 96 },
		{ 112 , 126 , 98 , 113 , -1 , 127 , -1 , 99 , -1 },
		{ 113 , 127 , 99 , 114 , 112 , 128 , 126 , 100 , 98 },
		{ 114 , 128 , 100 , 115 , 113 , 129 , 127 , 101 , 99 },
		{ 115 , 129 , 101 , 116 , 114 , 130 , 128 , 102 , 100 },
		{ 116 , 130 , 102 , 117 , 115 , 131 , 129 , 103 , 101 },
		{ 117 , 131 , 103 , 118 , 116 , 132 , 130 , 104 , 102 },
		{ 118 , 132 , 104 , 119 , 117 , 133 , 131 , 105 , 103 },
		{ 119 , 133 , 105 , 120 , 118 , 134 , 132 , 106 , 104 },
		{ 120 , 134 , 106 , 121 , 119 , 135 , 133 , 107 , 105 },
		{ 121 , 135 , 107 , 122 , 120 , 136 , 134 , 108 , 106 },
		{ 122 , 136 , 108 , 123 , 121 , 137 , 135 , 109 , 107 },
		{ 123 , 137 , 109 , 124 , 122 , 138 , 136 , 110 , 108 },
		{ 124 , 138 , 110 , 125 , 123 , 139 , 137 , 111 , 109 },
		{ 125 , 139 , 111 , -1 , 124 , -1 , 138 , -1 , 110 },
		{ 126 , 140 , 112 , 127 , -1 , 141 , -1 , 113 , -1 },
		{ 127 , 141 , 113 , 128 , 126 , 142 , 140 , 114 , 112 },
		{ 128 , 142 , 114 , 129 , 127 , 143 , 141 , 115 , 113 },
		{ 129 , 143 , 115 , 130 , 128 , 144 , 142 , 116 , 114 },
		{ 130 , 144 , 116 , 131 , 129 , 145 , 143 , 117 , 115 },
		{ 131 , 145 , 117 , 132 , 130 , 146 , 144 , 118 , 116 },
		{ 132 , 146 , 118 , 133 , 131 , 147 , 145 , 119 , 117 },
		{ 133 , 147 , 119 , 134 , 132 , 148 , 146 , 120 , 118 },
		{ 134 , 148 , 120 , 135 , 133 , 149 , 147 , 121 , 119 },
		{ 135 , 149 , 121 , 136 , 134 , 150 , 148 , 122 , 120 },
		{ 136 , 150 , 122 , 137 , 135 , 151 , 149 , 123 , 121 },
		{ 137 , 151 , 123 , 138 , 136 , 152 , 150 , 124 , 122 },
		{ 138 , 152 , 124 , 139 , 137 , 153 , 151 , 125 , 123 },
		{ 139 , 153 , 125 , -1 , 138 , -1 , 152 , -1 , 124 },
		{ 140 , 154 , 126 , 141 , -1 , 155 , -1 , 127 , -1 },
		{ 141 , 155 , 127 , 142 , 140 , 156 , 154 , 128 , 126 },
		{ 142 , 156 , 128 , 143 , 141 , 157 , 155 , 129 , 127 },
		{ 143 , 157 , 129 , 144 , 142 , 158 , 156 , 130 , 128 },
		{ 144 , 158 , 130 , 145 , 143 , 159 , 157 , 131 , 129 },
		{ 145 , 159 , 131 , 146 , 144 , 160 , 158 , 132 , 130 },
		{ 146 , 160 , 132 , 147 , 145 , 161 , 159 , 133 , 131 },
		{ 147 , 161 , 133 , 148 , 146 , 162 , 160 , 134 , 132 },
		{ 148 , 162 , 134 , 149 , 147 , 163 , 161 , 135 , 133 },
		{ 149 , 163 , 135 , 150 , 148 , 164 , 162 , 136 , 134 },
		{ 150 , 164 , 136 , 151 , 149 , 165 , 163 , 137 , 135 },
		{ 151 , 165 , 137 , 152 , 150 , 166 , 164 , 138 , 136 },
		{ 152 , 166 , 138 , 153 , 151 , 167 , 165 , 139 , 137 },
		{ 153 , 167 , 139 , -1 , 152 , -1 , 166 , -1 , 138 },
		{ 154 , 168 , 140 , 155 , -1 , 169 , -1 , 141 , -1 },
		{ 155 , 169 , 141 , 156 , 154 , 170 , 168 , 142 , 140 },
		{ 156 , 170 , 142 , 157 , 155 , 171 , 169 , 143 , 141 },
		{ 157 , 171 , 143 , 158 , 156 , 172 , 170 , 144 , 142 },
		{ 158 , 172 , 144 , 159 , 157 , 173 , 171 , 145 , 143 },
		{ 159 , 173 , 145 , 160 , 158 , 174 , 172 , 146 , 144 },
		{ 160 , 174 , 146 , 161 , 159 , 175 , 173 , 147 , 145 },
		{ 161 , 175 , 147 , 162 , 160 , 176 , 174 , 148 , 146 },
		{ 162 , 176 , 148 , 163 , 161 , 177 , 175 , 149 , 147 },
		{ 163 , 177 , 149 , 164 , 162 , 178 , 176 , 150 , 148 },
		{ 164 , 178 , 150 , 165 , 163 , 179 , 177 , 151 , 149 },
		{ 165 , 179 , 151 , 166 , 164 , 180 , 178 , 152 , 150 },
		{ 166 , 180 , 152 , 167 , 165 , 181 , 179 , 153 , 151 },
		{ 167 , 181 , 153 , -1 , 166 , -1 , 180 , -1 , 152 },
		{ 168 , 182 , 154 , 169 , -1 , 183 , -1 , 155 , -1 },
		{ 169 , 183 , 155 , 170 , 168 , 184 , 182 , 156 , 154 },
		{ 170 , 184 , 156 , 171 , 169 , 185 , 183 , 157 , 155 },
		{ 171 , 185 , 157 , 172 , 170 , 186 , 184 , 158 , 156 },
		{ 172 , 186 , 158 , 173 , 171 , 187 , 185 , 159 , 157 },
		{ 173 , 187 , 159 , 174 , 172 , 188 , 186 , 160 , 158 },
		{ 174 , 188 , 160 , 175 , 173 , 189 , 187 , 161 , 159 },
		{ 175 , 189 , 161 , 176 , 174 , 190 , 188 , 162 , 160 },
		{ 176 , 190 , 162 , 177 , 175 , 191 , 189 , 163 , 161 },
		{ 177 , 191 , 163 , 178 , 176 , 192 , 190 , 164 , 162 },
		{ 178 , 192 , 164 , 179 , 177 , 193 , 191 , 165 , 163 },
		{ 179 , 193 , 165 , 180 , 178 , 194 , 192 , 166 , 164 },
		{ 180 , 194 , 166 , 181 , 179 , 195 , 193 , 167 , 165 },
		{ 181 , 195 , 167 , -1 , 180 , -1 , 194 , -1 , 166 },
		{ 182 , 196 , 168 , 183 , -1 , 197 , -1 , 169 , -1 },
		{ 183 , 197 , 169 , 184 , 182 , 198 , 196 , 170 , 168 },
		{ 184 , 198 , 170 , 185 , 183 , 199 , 197 , 171 , 169 },
		{ 185 , 199 , 171 , 186 , 184 , 200 , 198 , 172 , 170 },
		{ 186 , 200 , 172 , 187 , 185 , 201 , 199 , 173 , 171 },
		{ 187 , 201 , 173 , 188 , 186 , 202 , 200 , 174 , 172 },
		{ 188 , 202 , 174 , 189 , 187 , 203 , 201 , 175 , 173 },
		{ 189 , 203 , 175 , 190 , 188 , 204 , 202 , 176 , 174 },
		{ 190 , 204 , 176 , 191 , 189 , 205 , 203 , 177 , 175 },
		{ 191 , 205 , 177 , 192 , 190 , 206 , 204 , 178 , 176 },
		{ 192 , 206 , 178 , 193 , 191 , 207 , 205 , 179 , 177 },
		{ 193 , 207 , 179 , 194 , 192 , 208 , 206 , 180 , 178 },
		{ 194 , 208 , 180 , 195 , 193 , 209 , 207 , 181 , 179 },
		{ 195 , 209 , 181 , -1 , 194 , -1 , 208 , -1 , 180 },
		{ 196 , 210 , 182 , 197 , -1 , 211 , -1 , 183 , -1 },
		{ 197 , 211 , 183 , 198 , 196 , 212 , 210 , 184 , 182 },
		{ 198 , 212 , 184 , 199 , 197 , 213 , 211 , 185 , 183 },
		{ 199 , 213 , 185 , 200 , 198 , 214 , 212 , 186 , 184 },
		{ 200 , 214 , 186 , 201 , 199 , 215 , 213 , 187 , 185 },
		{ 201 , 215 , 187 , 202 , 200 , 216 , 214 , 188 , 186 },
		{ 202 , 216 , 188 , 203 , 201 , 217 , 215 , 189 , 187 },
		{ 203 , 217 , 189 , 204 , 202 , 218 , 216 , 190 , 188 },
		{ 204 , 218 , 190 , 205 , 203 , 219 , 217 , 191 , 189 },
		{ 205 , 219 , 191 , 206 , 204 , 220 , 218 , 192 , 190 },
		{ 206 , 220 , 192 , 207 , 205 , 221 , 219 , 193 , 191 },
		{ 207 , 221 , 193 , 208 , 206 , 222 , 220 , 194 , 192 },
		{ 208 , 222 , 194 , 209 , 207 , 223 , 221 , 195 , 193 },
		{ 209 , 223 , 195 , -1 , 208 , -1 , 222 , -1 , 194 },
		{ 210 , 224 , 196 , 211 , -1 , 225 , -1 , 197 , -1 },
		{ 211 , 225 , 197 , 212 , 210 , 226 , 224 , 198 , 196 },
		{ 212 , 226 , 198 , 213 , 211 , 227 , 225 , 199 , 197 },
		{ 213 , 227 , 199 , 214 , 212 , 228 , 226 , 200 , 198 },
		{ 214 , 228 , 200 , 215 , 213 , 229 , 227 , 201 , 199 },
		{ 215 , 229 , 201 , 216 , 214 , 230 , 228 , 202 , 200 },
		{ 216 , 230 , 202 , 217 , 215 , 231 , 229 , 203 , 201 },
		{ 217 , 231 , 203 , 218 , 216 , 232 , 230 , 204 , 202 },
		{ 218 , 232 , 204 , 219 , 217 , 233 , 231 , 205 , 203 },
		{ 219 , 233 , 205 , 220 , 218 , 234 , 232 , 206 , 204 },
		{ 220 , 234 , 206 , 221 , 219 , 235 , 233 , 207 , 205 },
		{ 221 , 235 , 207 , 222 , 220 , 236 , 234 , 208 , 206 },
		{ 222 , 236 , 208 , 223 , 221 , 237 , 235 , 209 , 207 },
		{ 223 , 237 , 209 , -1 , 222 , -1 , 236 , -1 , 208 },
		{ 224 , 238 , 210 , 225 , -1 , 239 , -1 , 211 , -1 },
		{ 225 , 239 , 211 , 226 , 224 , 240 , 238 , 212 , 210 },
		{ 226 , 240 , 212 , 227 , 225 , 241 , 239 , 213 , 211 },
		{ 227 , 241 , 213 , 228 , 226 , 242 , 240 , 214 , 212 },
		{ 228 , 242 , 214 , 229 , 227 , 243 , 241 , 215 , 213 },
		{ 229 , 243 , 215 , 230 , 228 , 244 , 242 , 216 , 214 },
		{ 230 , 244 , 216 , 231 , 229 , 245 , 243 , 217 , 215 },
		{ 231 , 245 , 217 , 232 , 230 , 246 , 244 , 218 , 216 },
		{ 232 , 246 , 218 , 233 , 231 , 247 , 245 , 219 , 217 },
		{ 233 , 247 , 219 , 234 , 232 , 248 , 246 , 220 , 218 },
		{ 234 , 248 , 220 , 235 , 233 , 249 , 247 , 221 , 219 },
		{ 235 , 249 , 221 , 236 , 234 , 250 , 248 , 222 , 220 },
		{ 236 , 250 , 222 , 237 , 235 , 251 , 249 , 223 , 221 },
		{ 237 , 251 , 223 , -1 , 236 , -1 , 250 , -1 , 222 },
		{ 238 , 0 , 224 , 239 , -1 , 1 , -1 , 225 , -1 },
		{ 239 , 1 , 225 , 240 , 238 , 2 , 0 , 226 , 224 },
		{ 240 , 2 , 226 , 241 , 239 , 3 , 1 , 227 , 225 },
		{ 241 , 3 , 227 , 242 , 240 , 4 , 2 , 228 , 226 },
		{ 242 , 4 , 228 , 243 , 241 , 5 , 3 , 229 , 227 },
		{ 243 , 5 , 229 , 244 , 242 , 6 , 4 , 230 , 228 },
		{ 244 , 6 , 230 , 245 , 243 , 7 , 5 , 231 , 229 },
		{ 245 , 7 , 231 , 246 , 244 , 8 , 6 , 232 , 230 },
		{ 246 , 8 , 232 , 247 , 245 , 9 , 7 , 233 , 231 },
		{ 247 , 9 , 233 , 248 , 246 , 10 , 8 , 234 , 232 },
		{ 248 , 10 , 234 , 249 , 247 , 11 , 9 , 235 , 233 },
		{ 249 , 11 , 235 , 250 , 248 , 12 , 10 , 236 , 234 },
		{ 250 , 12 , 236 , 251 , 249 , 13 , 11 , 237 , 235 },
		{ 251 , 13 , 237 , -1 , 250 , -1 , 12 , -1 , 236 }
};

#endif