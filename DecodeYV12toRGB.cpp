#include "include\avisynth.h"
#include <windows.h>
#include <immintrin.h>


enum {
	CS_YUVA = 1 << 27,
	CS_BGR = 1 << 28,
	CS_YUV = 1 << 29,
	CS_INTERLEAVED = 1 << 30,
	CS_PLANAR = 1 << 31,

	CS_Shift_Sub_Width = 0,
	CS_Shift_Sub_Height = 8,
	CS_Shift_Sample_Bits = 16,

	CS_Sub_Width_Mask = 7 << CS_Shift_Sub_Width,
	CS_Sub_Width_1 = 3 << CS_Shift_Sub_Width, // YV24
	CS_Sub_Width_2 = 0 << CS_Shift_Sub_Width, // YV12, I420, YV16
	CS_Sub_Width_4 = 1 << CS_Shift_Sub_Width, // YUV9, YV411

	CS_VPlaneFirst = 1 << 3, // YV12, YV16, YV24, YV411, YUV9
	CS_UPlaneFirst = 1 << 4, // I420

	CS_Sub_Height_Mask = 7 << CS_Shift_Sub_Height,
	CS_Sub_Height_1 = 3 << CS_Shift_Sub_Height, // YV16, YV24, YV411
	CS_Sub_Height_2 = 0 << CS_Shift_Sub_Height, // YV12, I420
	CS_Sub_Height_4 = 1 << CS_Shift_Sub_Height, // YUV9

	CS_Sample_Bits_Mask = 7 << CS_Shift_Sample_Bits,
	CS_Sample_Bits_8 = 0 << CS_Shift_Sample_Bits,
	CS_Sample_Bits_10 = 5 << CS_Shift_Sample_Bits,
	CS_Sample_Bits_12 = 6 << CS_Shift_Sample_Bits,
	CS_Sample_Bits_14 = 7 << CS_Shift_Sample_Bits,
	CS_Sample_Bits_16 = 1 << CS_Shift_Sample_Bits,
	CS_Sample_Bits_32 = 2 << CS_Shift_Sample_Bits,

	CS_PLANAR_MASK = CS_PLANAR | CS_INTERLEAVED | CS_YUV | CS_BGR | CS_YUVA | CS_Sample_Bits_Mask
	| CS_Sub_Height_Mask | CS_Sub_Width_Mask,
	CS_PLANAR_FILTER = ~(CS_VPlaneFirst | CS_UPlaneFirst),

	CS_RGB_TYPE = 1 << 0,
	CS_RGBA_TYPE = 1 << 1,

	// Specific colorformats
	CS_UNKNOWN = 0,

	CS_BGR24 = CS_RGB_TYPE | CS_BGR | CS_INTERLEAVED,
	CS_BGR32 = CS_RGBA_TYPE | CS_BGR | CS_INTERLEAVED,
	CS_YUY2 = 1 << 2 | CS_YUV | CS_INTERLEAVED,
	//  CS_YV12  = 1<<3  Reserved
	//  CS_I420  = 1<<4  Reserved
	CS_RAW32 = 1 << 5 | CS_INTERLEAVED,

	//  YV12 must be 0xA000008 2.5 Baked API will see all new planar as YV12
	//  I420 must be 0xA000010

	CS_GENERIC_YUV420 = CS_PLANAR | CS_YUV | CS_VPlaneFirst | CS_Sub_Height_2 | CS_Sub_Width_2,  // 4:2:0 planar
	CS_GENERIC_YUV422 = CS_PLANAR | CS_YUV | CS_VPlaneFirst | CS_Sub_Height_1 | CS_Sub_Width_2,  // 4:2:2 planar
	CS_GENERIC_YUV444 = CS_PLANAR | CS_YUV | CS_VPlaneFirst | CS_Sub_Height_1 | CS_Sub_Width_1,  // 4:4:4 planar
	CS_GENERIC_Y = CS_PLANAR | CS_INTERLEAVED | CS_YUV,                                     // Y only (4:0:0)
	CS_GENERIC_RGBP = CS_PLANAR | CS_BGR | CS_RGB_TYPE,                                        // planar RGB. Though name is RGB but plane order G,B,R
	CS_GENERIC_RGBAP = CS_PLANAR | CS_BGR | CS_RGBA_TYPE,                                       // planar RGBA
	CS_GENERIC_YUVA420 = CS_PLANAR | CS_YUVA | CS_VPlaneFirst | CS_Sub_Height_2 | CS_Sub_Width_2, // 4:2:0:A planar
	CS_GENERIC_YUVA422 = CS_PLANAR | CS_YUVA | CS_VPlaneFirst | CS_Sub_Height_1 | CS_Sub_Width_2, // 4:2:2:A planar
	CS_GENERIC_YUVA444 = CS_PLANAR | CS_YUVA | CS_VPlaneFirst | CS_Sub_Height_1 | CS_Sub_Width_1, // 4:4:4:A planar

	CS_YV24 = CS_GENERIC_YUV444 | CS_Sample_Bits_8,  // YVU 4:4:4 planar
	CS_YV16 = CS_GENERIC_YUV422 | CS_Sample_Bits_8,  // YVU 4:2:2 planar
	CS_YV12 = CS_GENERIC_YUV420 | CS_Sample_Bits_8,  // YVU 4:2:0 planar
	CS_I420 = CS_PLANAR | CS_YUV | CS_Sample_Bits_8 | CS_UPlaneFirst | CS_Sub_Height_2 | CS_Sub_Width_2,  // YUV 4:2:0 planar
	CS_IYUV = CS_I420,
	CS_YUV9 = CS_PLANAR | CS_YUV | CS_Sample_Bits_8 | CS_VPlaneFirst | CS_Sub_Height_4 | CS_Sub_Width_4,  // YUV 4:1:0 planar
	CS_YV411 = CS_PLANAR | CS_YUV | CS_Sample_Bits_8 | CS_VPlaneFirst | CS_Sub_Height_1 | CS_Sub_Width_4,  // YUV 4:1:1 planar

	CS_Y8 = CS_GENERIC_Y | CS_Sample_Bits_8,                                                            // Y   4:0:0 planar

	//-------------------------
	// AVS16: new planar constants go live! Experimental PF 160613
	// 10-12-14 bit + planar RGB + BRG48/64 160725

	CS_YUV444P10 = CS_GENERIC_YUV444 | CS_Sample_Bits_10, // YUV 4:4:4 10bit samples
	CS_YUV422P10 = CS_GENERIC_YUV422 | CS_Sample_Bits_10, // YUV 4:2:2 10bit samples
	CS_YUV420P10 = CS_GENERIC_YUV420 | CS_Sample_Bits_10, // YUV 4:2:0 10bit samples
	CS_Y10 = CS_GENERIC_Y | CS_Sample_Bits_10,            // Y   4:0:0 10bit samples

	CS_YUV444P12 = CS_GENERIC_YUV444 | CS_Sample_Bits_12, // YUV 4:4:4 12bit samples
	CS_YUV422P12 = CS_GENERIC_YUV422 | CS_Sample_Bits_12, // YUV 4:2:2 12bit samples
	CS_YUV420P12 = CS_GENERIC_YUV420 | CS_Sample_Bits_12, // YUV 4:2:0 12bit samples
	CS_Y12 = CS_GENERIC_Y | CS_Sample_Bits_12,            // Y   4:0:0 12bit samples

	CS_YUV444P14 = CS_GENERIC_YUV444 | CS_Sample_Bits_14, // YUV 4:4:4 14bit samples
	CS_YUV422P14 = CS_GENERIC_YUV422 | CS_Sample_Bits_14, // YUV 4:2:2 14bit samples
	CS_YUV420P14 = CS_GENERIC_YUV420 | CS_Sample_Bits_14, // YUV 4:2:0 14bit samples
	CS_Y14 = CS_GENERIC_Y | CS_Sample_Bits_14,            // Y   4:0:0 14bit samples

	CS_YUV444P16 = CS_GENERIC_YUV444 | CS_Sample_Bits_16, // YUV 4:4:4 16bit samples
	CS_YUV422P16 = CS_GENERIC_YUV422 | CS_Sample_Bits_16, // YUV 4:2:2 16bit samples
	CS_YUV420P16 = CS_GENERIC_YUV420 | CS_Sample_Bits_16, // YUV 4:2:0 16bit samples
	CS_Y16 = CS_GENERIC_Y | CS_Sample_Bits_16,            // Y   4:0:0 16bit samples

	// 32 bit samples (float)
	CS_YUV444PS = CS_GENERIC_YUV444 | CS_Sample_Bits_32,  // YUV 4:4:4 32bit samples
	CS_YUV422PS = CS_GENERIC_YUV422 | CS_Sample_Bits_32,  // YUV 4:2:2 32bit samples
	CS_YUV420PS = CS_GENERIC_YUV420 | CS_Sample_Bits_32,  // YUV 4:2:0 32bit samples
	CS_Y32 = CS_GENERIC_Y | CS_Sample_Bits_32,            // Y   4:0:0 32bit samples

	// RGB packed
	CS_BGR48 = CS_RGB_TYPE | CS_BGR | CS_INTERLEAVED | CS_Sample_Bits_16, // BGR 3x16 bit
	CS_BGR64 = CS_RGBA_TYPE | CS_BGR | CS_INTERLEAVED | CS_Sample_Bits_16, // BGR 4x16 bit
	// no packed 32 bit (float) support for these legacy types

	// RGB planar
	CS_RGBP = CS_GENERIC_RGBP | CS_Sample_Bits_8,  // Planar RGB 8 bit samples
	CS_RGBP8 = CS_GENERIC_RGBP | CS_Sample_Bits_8,  // Planar RGB 8 bit samples
	CS_RGBP10 = CS_GENERIC_RGBP | CS_Sample_Bits_10, // Planar RGB 10bit samples
	CS_RGBP12 = CS_GENERIC_RGBP | CS_Sample_Bits_12, // Planar RGB 12bit samples
	CS_RGBP14 = CS_GENERIC_RGBP | CS_Sample_Bits_14, // Planar RGB 14bit samples
	CS_RGBP16 = CS_GENERIC_RGBP | CS_Sample_Bits_16, // Planar RGB 16bit samples
	CS_RGBPS = CS_GENERIC_RGBP | CS_Sample_Bits_32, // Planar RGB 32bit samples

	// RGBA planar
	CS_RGBAP = CS_GENERIC_RGBAP | CS_Sample_Bits_8,  // Planar RGBA 8 bit samples
	CS_RGBAP8 = CS_GENERIC_RGBAP | CS_Sample_Bits_8,  // Planar RGBA 8 bit samples
	CS_RGBAP10 = CS_GENERIC_RGBAP | CS_Sample_Bits_10, // Planar RGBA 10bit samples
	CS_RGBAP12 = CS_GENERIC_RGBAP | CS_Sample_Bits_12, // Planar RGBA 12bit samples
	CS_RGBAP14 = CS_GENERIC_RGBAP | CS_Sample_Bits_14, // Planar RGBA 14bit samples
	CS_RGBAP16 = CS_GENERIC_RGBAP | CS_Sample_Bits_16, // Planar RGBA 16bit samples
	CS_RGBAPS = CS_GENERIC_RGBAP | CS_Sample_Bits_32, // Planar RGBA 32bit samples

	// Planar YUVA
	CS_YUVA444 = CS_GENERIC_YUVA444 | CS_Sample_Bits_8,  // YUVA 4:4:4 8bit samples
	CS_YUVA422 = CS_GENERIC_YUVA422 | CS_Sample_Bits_8,  // YUVA 4:2:2 8bit samples
	CS_YUVA420 = CS_GENERIC_YUVA420 | CS_Sample_Bits_8,  // YUVA 4:2:0 8bit samples

	CS_YUVA444P10 = CS_GENERIC_YUVA444 | CS_Sample_Bits_10, // YUVA 4:4:4 10bit samples
	CS_YUVA422P10 = CS_GENERIC_YUVA422 | CS_Sample_Bits_10, // YUVA 4:2:2 10bit samples
	CS_YUVA420P10 = CS_GENERIC_YUVA420 | CS_Sample_Bits_10, // YUVA 4:2:0 10bit samples

	CS_YUVA444P12 = CS_GENERIC_YUVA444 | CS_Sample_Bits_12, // YUVA 4:4:4 12bit samples
	CS_YUVA422P12 = CS_GENERIC_YUVA422 | CS_Sample_Bits_12, // YUVA 4:2:2 12bit samples
	CS_YUVA420P12 = CS_GENERIC_YUVA420 | CS_Sample_Bits_12, // YUVA 4:2:0 12bit samples

	CS_YUVA444P14 = CS_GENERIC_YUVA444 | CS_Sample_Bits_14, // YUVA 4:4:4 14bit samples
	CS_YUVA422P14 = CS_GENERIC_YUVA422 | CS_Sample_Bits_14, // YUVA 4:2:2 14bit samples
	CS_YUVA420P14 = CS_GENERIC_YUVA420 | CS_Sample_Bits_14, // YUVA 4:2:0 14bit samples

	CS_YUVA444P16 = CS_GENERIC_YUVA444 | CS_Sample_Bits_16, // YUVA 4:4:4 16bit samples
	CS_YUVA422P16 = CS_GENERIC_YUVA422 | CS_Sample_Bits_16, // YUVA 4:2:2 16bit samples
	CS_YUVA420P16 = CS_GENERIC_YUVA420 | CS_Sample_Bits_16, // YUVA 4:2:0 16bit samples

	CS_YUVA444PS = CS_GENERIC_YUVA444 | CS_Sample_Bits_32,  // YUVA 4:4:4 32bit samples
	CS_YUVA422PS = CS_GENERIC_YUVA422 | CS_Sample_Bits_32,  // YUVA 4:2:2 32bit samples
	CS_YUVA420PS = CS_GENERIC_YUVA420 | CS_Sample_Bits_32,  // YUVA 4:2:0 32bit samples

};


void Convert(PVideoFrame dst, PVideoFrame src, VideoInfo vi_dst, VideoInfo vi_src, int threads, int cpuFlags)
{
	auto srcp_Y = src->GetReadPtr(PLANAR_Y);
	auto srcp_U = src->GetReadPtr(PLANAR_U);
	auto srcp_V = src->GetReadPtr(PLANAR_V);


	auto height = src->GetHeight(PLANAR_Y);
	auto row_size_Y = src->GetRowSize(PLANAR_Y);
	auto row_size_U = src->GetRowSize(PLANAR_U);
	auto row_size_V = src->GetRowSize(PLANAR_V);
	auto src_pitch_Y = src->GetPitch(PLANAR_Y);
	auto src_pitch_U = src->GetPitch(PLANAR_U);
	auto src_pitch_V = src->GetPitch(PLANAR_V);

	/*
	auto dstp_R = dst->GetWritePtr(PLANAR_R_ALIGNED);
	auto dstp_G = dst->GetWritePtr(PLANAR_G_ALIGNED);
	auto dstp_B = dst->GetWritePtr(PLANAR_B_ALIGNED);
	*/
	auto dstp_R = dst->GetWritePtr(4);
	auto dstp_G = dst->GetWritePtr(6);
	auto dstp_B = dst->GetWritePtr(2);


	auto dst_pitch_R = dst->GetPitch(PLANAR_R);
	auto dst_pitch_G = dst->GetPitch(PLANAR_G);
	auto dst_pitch_B = dst->GetPitch(PLANAR_B);


#pragma omp parallel for num_threads(threads)
	for (int y = 0; y < height; y++)
	{
		unsigned char* l_dstp_R = dstp_R + y * dst_pitch_R;
		unsigned char* l_dstp_G = dstp_G + y * dst_pitch_G;
		unsigned char* l_dstp_B = dstp_B + y * dst_pitch_B;

		unsigned char* l_srcp_Y = (unsigned char*)srcp_Y + y * src_pitch_Y;
		unsigned char* l_srcp_U = (unsigned char*)srcp_U + (y >> 1) * src_pitch_U;
		unsigned char* l_srcp_V = (unsigned char*)srcp_V + (y >> 1) * src_pitch_V;


			if (cpuFlags & CPUF_AVX2) // use AVX2
			{
				const int col64 = row_size_Y - (row_size_Y % 64); // load 64 Y samples and 

				__m256i ymm_w128 = _mm256_set1_epi16(128);

				__m256i ymm_w73 = _mm256_set1_epi16(73); // Kr of 601 ? , div64
				__m256i ymm_w130 = _mm256_set1_epi16(130); // Kb of 601 ? , div64

				__m256i ymm_w25 = _mm256_set1_epi16(25); // Kgu of 601 ? , div64
				__m256i ymm_w37 = _mm256_set1_epi16(37); // Kgv of 601 ? , div64

				for (int col = 0; col < col64; col += 64)
				{
					__m256i ymm0_Y0 = _mm256_lddqu_si256((const __m256i*)l_srcp_Y); // better align start addr with pre-conversion of 32-bytes aligned (if exist) and use load_ps
					__m256i ymm1_Y1 = _mm256_lddqu_si256((const __m256i*)(l_srcp_Y + 32));
					__m256i ymm2_U = _mm256_lddqu_si256((const __m256i*)(l_srcp_U));
					__m256i ymm3_V = _mm256_lddqu_si256((const __m256i*)(l_srcp_V));

					__m256i ymm_Y0_16l = _mm256_unpacklo_epi8(ymm0_Y0, _mm256_setzero_si256());
					__m256i ymm_Y1_16l = _mm256_unpacklo_epi8(ymm1_Y1, _mm256_setzero_si256());

					__m256i ymm_Y0_16h = _mm256_unpackhi_epi8(ymm0_Y0, _mm256_setzero_si256());
					__m256i ymm_Y1_16h = _mm256_unpackhi_epi8(ymm1_Y1, _mm256_setzero_si256());

					__m256i ymm_U_dl = _mm256_unpacklo_epi8(ymm2_U, ymm2_U);
					__m256i ymm_V_dl = _mm256_unpacklo_epi8(ymm3_V, ymm3_V);

					__m256i ymm_U_dh = _mm256_unpackhi_epi8(ymm2_U, ymm2_U);
					__m256i ymm_V_dh = _mm256_unpackhi_epi8(ymm3_V, ymm3_V);

					__m256i ymm_U_dl16l = _mm256_unpacklo_epi8(ymm_U_dl, _mm256_setzero_si256());
					__m256i ymm_V_dl16l = _mm256_unpacklo_epi8(ymm_V_dl, _mm256_setzero_si256());

					__m256i ymm_U_dl16h = _mm256_unpackhi_epi8(ymm_U_dl, _mm256_setzero_si256());
					__m256i ymm_V_dl16h = _mm256_unpackhi_epi8(ymm_V_dl, _mm256_setzero_si256());

					__m256i ymm_U_dh16l = _mm256_unpacklo_epi8(ymm_U_dh, _mm256_setzero_si256());
					__m256i ymm_V_dh16l = _mm256_unpacklo_epi8(ymm_V_dh, _mm256_setzero_si256());

					__m256i ymm_U_dh16h = _mm256_unpackhi_epi8(ymm_U_dh, _mm256_setzero_si256());
					__m256i ymm_V_dh16h = _mm256_unpackhi_epi8(ymm_V_dh, _mm256_setzero_si256());


					ymm_U_dl16l = _mm256_sub_epi16(ymm_U_dl16l, ymm_w128); // superscalarity of 3 at IceLake/SkyLake
					ymm_V_dl16l = _mm256_sub_epi16(ymm_V_dl16l, ymm_w128);

					ymm_U_dl16h = _mm256_sub_epi16(ymm_U_dl16h, ymm_w128);
					ymm_V_dl16h = _mm256_sub_epi16(ymm_V_dl16h, ymm_w128);

					ymm_U_dh16l = _mm256_sub_epi16(ymm_U_dh16l, ymm_w128);
					ymm_V_dh16l = _mm256_sub_epi16(ymm_V_dh16l, ymm_w128);

					ymm_U_dh16h = _mm256_sub_epi16(ymm_U_dh16h, ymm_w128);
					ymm_V_dh16h = _mm256_sub_epi16(ymm_V_dh16h, ymm_w128);


					__m256i ymm_V_dl16l_addR = _mm256_mullo_epi16(ymm_V_dl16l, ymm_w73);
					__m256i ymm_V_dl16h_addR = _mm256_mullo_epi16(ymm_V_dl16h, ymm_w73);

					__m256i ymm_V_dh16l_addR = _mm256_mullo_epi16(ymm_V_dh16l, ymm_w73);
					__m256i ymm_V_dh16h_addR = _mm256_mullo_epi16(ymm_V_dh16h, ymm_w73);

					
					__m256i ymm_U_dl16l_addB = _mm256_mullo_epi16(ymm_U_dl16l, ymm_w130);
					__m256i ymm_U_dl16h_addB = _mm256_mullo_epi16(ymm_U_dl16h, ymm_w130);

					__m256i ymm_U_dh16l_addB = _mm256_mullo_epi16(ymm_U_dh16l, ymm_w130);
					__m256i ymm_U_dh16h_addB = _mm256_mullo_epi16(ymm_U_dh16h, ymm_w130);

					ymm_V_dl16l_addR = _mm256_srai_epi16(ymm_V_dl16l_addR, 6);
					ymm_V_dl16h_addR = _mm256_srai_epi16(ymm_V_dl16h_addR, 6);

					ymm_V_dh16l_addR = _mm256_srai_epi16(ymm_V_dh16l_addR, 6);
					ymm_V_dh16h_addR = _mm256_srai_epi16(ymm_V_dh16h_addR, 6);


					ymm_U_dl16l_addB = _mm256_srai_epi16(ymm_U_dl16l_addB, 6);
					ymm_U_dl16h_addB = _mm256_srai_epi16(ymm_U_dl16h_addB, 6);

					ymm_U_dh16l_addB = _mm256_srai_epi16(ymm_U_dh16l_addB, 6);
					ymm_U_dh16h_addB = _mm256_srai_epi16(ymm_U_dh16h_addB, 6);


					// R
					__m256i ymm_R_0_16l = _mm256_add_epi16(ymm_Y0_16l, ymm_V_dl16l_addR);
					__m256i ymm_R_0_16h = _mm256_add_epi16(ymm_Y0_16h, ymm_V_dl16h_addR);

					__m256i ymm_R_1_16l = _mm256_add_epi16(ymm_Y1_16l, ymm_V_dh16l_addR);
					__m256i ymm_R_1_16h = _mm256_add_epi16(ymm_Y1_16h, ymm_V_dh16h_addR);

					// B
					__m256i ymm_B_0_16l = _mm256_add_epi16(ymm_Y0_16l, ymm_U_dl16l_addB);
					__m256i ymm_B_0_16h = _mm256_add_epi16(ymm_Y0_16h, ymm_U_dl16h_addB);

					__m256i ymm_B_1_16l = _mm256_add_epi16(ymm_Y1_16l, ymm_U_dh16l_addB);
					__m256i ymm_B_1_16h = _mm256_add_epi16(ymm_Y1_16h, ymm_U_dh16h_addB);

					//G
					__m256i ymm_V_dl16l_subG = _mm256_mullo_epi16(ymm_V_dl16l, ymm_w37);
					__m256i ymm_V_dl16h_subG = _mm256_mullo_epi16(ymm_V_dl16h, ymm_w37);

					__m256i ymm_V_dh16l_subG = _mm256_mullo_epi16(ymm_V_dh16l, ymm_w37);
					__m256i ymm_V_dh16h_subG = _mm256_mullo_epi16(ymm_V_dh16h, ymm_w37);


					__m256i ymm_U_dl16l_subG = _mm256_mullo_epi16(ymm_U_dl16l, ymm_w25);
					__m256i ymm_U_dl16h_subG = _mm256_mullo_epi16(ymm_U_dl16h, ymm_w25);

					__m256i ymm_U_dh16l_subG = _mm256_mullo_epi16(ymm_U_dh16l, ymm_w25);
					__m256i ymm_U_dh16h_subG = _mm256_mullo_epi16(ymm_U_dh16h, ymm_w25);


					ymm_V_dl16l_subG = _mm256_srai_epi16(ymm_V_dl16l_subG, 6);
					ymm_V_dl16h_subG = _mm256_srai_epi16(ymm_V_dl16h_subG, 6);

					ymm_V_dh16l_subG = _mm256_srai_epi16(ymm_V_dh16l_subG, 6);
					ymm_V_dh16h_subG = _mm256_srai_epi16(ymm_V_dh16h_subG, 6);


					ymm_U_dl16l_subG = _mm256_srai_epi16(ymm_U_dl16l_subG, 6);
					ymm_U_dl16h_subG = _mm256_srai_epi16(ymm_U_dl16h_subG, 6);

					ymm_U_dh16l_subG = _mm256_srai_epi16(ymm_U_dh16l_subG, 6);
					ymm_U_dh16h_subG = _mm256_srai_epi16(ymm_U_dh16h_subG, 6);


					__m256i ymm_G_0_16l = _mm256_sub_epi16(ymm_Y0_16l, ymm_V_dl16l_subG);
					__m256i ymm_G_0_16h = _mm256_sub_epi16(ymm_Y0_16h, ymm_V_dl16h_subG);

					ymm_G_0_16l = _mm256_sub_epi16(ymm_G_0_16l, ymm_U_dl16l_subG);
					ymm_G_0_16h = _mm256_sub_epi16(ymm_G_0_16h, ymm_U_dl16h_subG);

					__m256i ymm_G_1_16l = _mm256_sub_epi16(ymm_Y1_16l, ymm_V_dh16l_subG);
					__m256i ymm_G_1_16h = _mm256_sub_epi16(ymm_Y1_16h, ymm_V_dh16h_subG);

					ymm_G_1_16l = _mm256_sub_epi16(ymm_G_0_16l, ymm_U_dh16l_subG);
					ymm_G_1_16h = _mm256_sub_epi16(ymm_G_0_16h, ymm_U_dh16h_subG);

					//pack 16bit to 8bit
					__m256i ymm_R_0_8 = _mm256_packus_epi16(ymm_R_0_16l, ymm_R_0_16h);
					__m256i ymm_R_1_8 = _mm256_packus_epi16(ymm_R_1_16l, ymm_R_1_16h);

					__m256i ymm_G_0_8 = _mm256_packus_epi16(ymm_G_0_16l, ymm_G_0_16h);
					__m256i ymm_G_1_8 = _mm256_packus_epi16(ymm_G_1_16l, ymm_G_1_16h);

					__m256i ymm_B_0_8 = _mm256_packus_epi16(ymm_B_0_16l, ymm_B_0_16h);
					__m256i ymm_B_1_8 = _mm256_packus_epi16(ymm_B_1_16l, ymm_B_1_16h);

					_mm256_store_si256((__m256i*)l_dstp_R, ymm_R_0_8);
					_mm256_store_si256((__m256i*)(l_dstp_R + 32), ymm_R_1_8);

					_mm256_store_si256((__m256i*)l_dstp_G, ymm_G_0_8);
					_mm256_store_si256((__m256i*)(l_dstp_G + 32), ymm_G_1_8);

					_mm256_store_si256((__m256i*)l_dstp_B, ymm_B_0_8);
					_mm256_store_si256((__m256i*)(l_dstp_B + 32), ymm_B_1_8);


					l_srcp_Y += 64; // in bytes
					l_srcp_U += 32;
					l_srcp_V += 32;

					l_dstp_R += 64;
					l_dstp_G += 64;
					l_dstp_B += 64;
				}

				// last cols
				for (int col = col64; col < row_size_Y; ++col)
				{
				}
			}
	}
}

/*
template void Invert<uint8_t>(unsigned char* _srcp, unsigned char* _dstp, int src_pitch, int dst_pitch, int height, int row_size, int bits, int threads, int cpuFlags);
template void Invert<uint16_t>(unsigned char* _srcp, unsigned char* _dstp, int src_pitch, int dst_pitch, int height, int row_size, int bits, int threads, int cpuFlags);
template void Invert<float>(unsigned char* _srcp, unsigned char* _dstp, int src_pitch, int dst_pitch, int height, int row_size, int bits, int threads, int cpuFlags);
*/

class DecodeYV12toRGB : public GenericVideoFilter
{
	int threads;
	int _cpuFlags;

public:
	DecodeYV12toRGB(PClip _child, int threads_, IScriptEnvironment* env) : GenericVideoFilter(_child), threads(threads_)
	{
		_cpuFlags = env->GetCPUFlags();
		vi.pixel_type = CS_RGBP8;//CS_BGR24;
	}

	PVideoFrame __stdcall GetFrame(int n, IScriptEnvironment* env)
	{
		PVideoFrame dst = env->NewVideoFrame(vi);

		VideoInfo vi_src = child->GetVideoInfo();

		PVideoFrame src = child->GetFrame(n, env);

		if (vi_src.ComponentSize() == 1)
		{
			Convert(dst, src, vi, vi_src,  threads, _cpuFlags);
		}
		else
			env->ThrowError("DecodeYV12toRGB: Only 8bit input supported.");
		return dst;
	}
};


AVSValue __cdecl Create_Decode(AVSValue args, void* user_data, IScriptEnvironment* env)
{
	return new DecodeYV12toRGB(args[0].AsClip(), args[1].AsInt(1), env);
}

const AVS_Linkage* AVS_linkage = 0;

extern "C" __declspec(dllexport) const char* __stdcall AvisynthPluginInit3(IScriptEnvironment * env, const AVS_Linkage* const vectors)
{
	AVS_linkage = vectors;
	env->AddFunction("DecodeYV12toRGB", "c[threads]i", Create_Decode, 0);

	return "Decode YV12 to RGB sample plugin";
}


