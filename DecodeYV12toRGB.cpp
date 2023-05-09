#include "include\avisynth.h"
#include <windows.h>
#include <immintrin.h>


enum {
	CS_BGR = 1 << 28,
	CS_INTERLEAVED = 1 << 30,

	CS_RGB_TYPE = 1 << 0,
	CS_RGBA_TYPE = 1 << 1,

	CS_BGR24 = CS_RGB_TYPE | CS_BGR | CS_INTERLEAVED

};

//void Invert(unsigned char* _srcp, unsigned char* _dstp, int src_pitch, int dst_pitch, int height, int row_size, int bits, int threads, int cpuFlags)
void Invert(PVideoFrame dst, PVideoFrame src, VideoInfo vi_dst, VideoInfo vi_src, int threads, int cpuFlags)
{
	auto srcp_Y = src->GetReadPtr(PLANAR_Y);
	auto srcp_U = src->GetReadPtr(PLANAR_U);
	auto srcp_V = src->GetReadPtr(PLANAR_V);

	auto dstp = dst->GetWritePtr();
	auto height = src->GetHeight(PLANAR_Y);
	auto row_size_Y = src->GetRowSize(PLANAR_Y) / vi_src.ComponentSize();
	auto row_size_U = src->GetRowSize(PLANAR_U) / vi_src.ComponentSize();
	auto row_size_V = src->GetRowSize(PLANAR_V) / vi_src.ComponentSize();
	auto src_pitch_Y = src->GetPitch(PLANAR_Y) / vi_src.ComponentSize();
	auto src_pitch_U = src->GetPitch(PLANAR_U) / vi_src.ComponentSize();
	auto src_pitch_V = src->GetPitch(PLANAR_V) / vi_src.ComponentSize();
	auto dst_pitch = dst->GetPitch() / vi_dst.ComponentSize();


#pragma omp parallel for num_threads(threads)
	for (int y = 0; y < height; y++)
	{
		unsigned char* l_dstp = dstp + y * dst_pitch;
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
					__m256i ymm_R_0_8 = _mm256_packs_epi16(ymm_R_0_16l, ymm_R_0_16h);
					__m256i ymm_R_1_8 = _mm256_packs_epi16(ymm_R_1_16l, ymm_R_1_16h);

					__m256i ymm_G_0_8 = _mm256_packs_epi16(ymm_G_0_16l, ymm_G_0_16h);
					__m256i ymm_G_1_8 = _mm256_packs_epi16(ymm_G_1_16l, ymm_G_1_16h);

					__m256i ymm_B_0_8 = _mm256_packs_epi16(ymm_B_0_16l, ymm_B_0_16h);
					__m256i ymm_B_1_8 = _mm256_packs_epi16(ymm_B_1_16l, ymm_B_1_16h);


					l_srcp_Y += 64; // in bytes
					l_srcp_U += 32;
					l_srcp_V += 32;
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
		vi.pixel_type = CS_BGR24;
	}

	PVideoFrame __stdcall GetFrame(int n, IScriptEnvironment* env)
	{
		PVideoFrame dst = env->NewVideoFrame(vi);

		VideoInfo vi_src = child->GetVideoInfo();

		PVideoFrame src = child->GetFrame(n, env);

		/*
		auto srcp = src->GetReadPtr(PLANAR_Y);
		auto dstp = dst->GetWritePtr(PLANAR_Y);
		auto height = src->GetHeight(PLANAR_Y);
		auto row_size = src->GetRowSize(PLANAR_Y) / vi_src.ComponentSize();
		auto src_pitch = src->GetPitch(PLANAR_Y) / vi_src.ComponentSize();
		auto dst_pitch = dst->GetPitch(PLANAR_Y) / vi.ComponentSize();
		*/

		if (vi_src.ComponentSize() == 1)
		{
//			Invert((uint8_t*)srcp, dstp, src_pitch, dst_pitch, height, row_size, vi_src.BitsPerComponent(), threads, _cpuFlags);
			Invert(dst, src, vi, vi_src,  threads, _cpuFlags);
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
	env->AddFunction("DecodeYV12toRGB24", "c[threads]i", Create_Decode, 0);

	return "Decode YV12 to RGB sample plugin";
}


