#pragma once
#include <stdint.h>
#ifdef __cplusplus
extern "C" {
#endif
	typedef void* BeamformingGPUHandle;
	//只做波束合成的初始化
	__declspec(dllexport) BeamformingGPUHandle initializeBeamformingGPU(
		int Demod_AFE_Dynamic,
		int rx_num,
		int Channel,
		int RxChannel,
		int element_num,
		float element_pitch,
		int sample_num,
		int bf_sample_num,
		int BeamN,
		int ValidBeamN,
		float sos,
		float fs,
		float B_fs,
		float focus_depth,
		float Js_Min_Aper,
		float Js_Max_Aper,
		float Js_Max_Fn,
		float Js_Min_Fn,
		float Js_Fn_depth,
		float Js_fn_setp,
		float Js_Fn_value,
		float* allBeamN, int allBeamN_size,
		float* allValidBeamN, int allValidBeamN_size,
		float* allBeamNIndex, int allBeamNIndex_size,
		float* allbfframeIndex, int allbfframeIndex_size,
		float* allBFtype, int allBFtype_size,
		float* allfsfocus, int allfsfocus_size,
		float* allfslinex, int allfslinex_size,
		float* allfslinez, int allfslinez_size,
		float* Steering, int Steering_size,
		float* plae_wave_steer, int plae_wave_steer_size,
		float* beamx, int beamx_size,
		float* beamz, int beamz_size,
		float* cstartoffset, int cstartoffset_size,
		float* Jswin, int Jswin_size,
		float* udm_L,int udm_L_size,
		float* fs_focus,
		int fs_focus_size,
		float* channel_map, int channel_map_size,
		float* element_pos_xz, int element_pos_xz_size,

		float Demod_MaxPoint,
		float blackholethre,
		float ImageStart,
		float UI_depth,
		float* Demodulate_FreDepth,
		float* Demodulate_FreValue,
		float* Demod_FilterDepth,
		float* Demod_FilterCutoff,
		const char* Demodulate_FilterType,
		int Demod_OrderFactor,
		float wininfo_alpha
	);

	__declspec(dllexport) BeamformingGPUHandle initializespectralGPU(
		int prf,
		int numperfile,
		int buffer_num,
		int fft_period,
		int lag,
		float t_axis_span,
		int t_idx,
		int t_buffer_idx,
		const char* probe_type,
		int Demod_AFE_Dynamic,
		int rx_num,
		int Channel,
		int RxChannel,
		int element_num,
		float element_pitch,
		int adc_head,
		int bf_sample_num,
		int BeamN,
		float sos,
		float fs,
		float fc, //发射频率
		float focus_depth,
		float Js_Min_Aper,
		float Js_Max_Aper,
		float Js_Max_Fn,
		float Js_Min_Fn,
		float Js_Fn_depth,
		float Js_fn_setp,
		float Js_Fn_value,
		float Fsweep_step,
		float Fsweep_start,
		float Fsweep_end,
		float Roi_start,
		float* plae_wave_steer, int plae_wave_steer_size,
		float* beamx, int beamx_size,
		float* beamz, int beamz_size,
		float* cstartoffset, int cstartoffset_size,
		float* Jswin, int Jswin_size,
		float* udm_L, int udm_L_size,
		float* fs_focus, int fs_focus_size,
		float* channel_map, int channel_map_size,
		float* element_pos_xz, int element_pos_xz_size,
		float* filter_coef, int filter_size,
		float* rx_apod, int rx_apod_size,
		float* rx_delay, int rx_delay_size,
		float* tx_delay, int tx_delay_size,
		float* x_loc_all, float* z_loc_all, int loc_num,
		float* iir_b, float* iir_a, float* iir_zi, int iir_order
	);
	 
	__declspec(dllexport) 	BeamformingGPUHandle  initializeAllGPU(
		const char* probe_type,
		int Demod_AFE_Dynamic,
		int rx_num,
		int Channel,
		int RxChannel,
		int element_num,
		float element_pitch,
		int adc_head,
		int bf_sample_num,
		int BeamN,
		float sos,
		float fs,
		float focus_depth,
		float Js_Min_Aper,
		float Js_Max_Aper,
		float Js_Max_Fn,
		float Js_Min_Fn,
		float Js_Fn_depth,
		float Js_fn_setp,
		float Js_Fn_value,
		float Fsweep_step,
		float Fsweep_start,
		float Fsweep_end,
		float Roi_start,
		float* plae_wave_steer, int plae_wave_steer_size,
		float* beamx, int beamx_size,
		float* beamz, int beamz_size,
		float* cstartoffset, int cstartoffset_size,
		float* Jswin, int Jswin_size,
		float* udm_L, int udm_L_size,
		float* fs_focus, int fs_focus_size,
		float* channel_map, int channel_map_size,
		float* element_pos_xz, int element_pos_xz_size,
		int ValidBeamN,
		int rx_line_size,
		float* rx_line,
		float SmoothCoef,
		float UI_FrameRate, 
		float B_fs,
		int Demod_MaxPoint,
		float ImageStart,
		float UI_depth,
		int dr_input,
		int dr_output,
		int dr_drange,
		int axialGain_size,
		float* axialGain,
		float xuanniuGain,
		int spa_smoothcoef_size,
		float* spa_smoothcoef,
		float* Demodulate_FreDepth,
		float* Demodulate_FreValue,
		const char* Demodulate_FilterType,
		float blackholethre,
		float decifs,
		float* Demod_FilterDepth,
		float* Demod_FilterCutoff,
		int  Demod_OrderFactor,
		float wininfo_alpha
	);

	__declspec(dllexport) BeamformingGPUHandle initializeULMGPU(
		int numberOfParticles,
		int res,
		int SVD_cutoff_min,		int SVD_cutoff_max,
		int max_linking_distance,
		int min_length,
		int fwhm_x,		int fwhm_z,
		int max_gap_closing,
		int Nz,		int Nx,		int NbFrames,
		float scalez,		float scalex,		float scalet,
		float CutoffFreq_min,		float CutoffFreq_max,
		float framerate,
		bool flag_vel_norm, bool flag_vel_z
	);
	
	__declspec(dllexport) 	BeamformingGPUHandle  initializeAllGPU_Test();
	//实时调用
	//return: <0  失败 >0 表示没有形成一张图,比如平面波的前面几帧  == 0 成功
	__declspec(dllexport) int processDataBeamformingGPU(BeamformingGPUHandle handle, int8_t* input, int size, float* output);
	// 中处理接口,输入:1.RF数据 2.线密度 3.模式 4.中处理数据  返回:,图像宽度
	__declspec(dllexport) int processDataMidProcessGPU(BeamformingGPUHandle handle, float* input,int noLine,int modetype, float* output);

	__declspec(dllexport) int processSpectralDataBeamformingandPostGPU(BeamformingGPUHandle handle, int8_t* input, 
		int size, float* output_b, float* output_sd, int bag_idx);
	
	//后处理
	__declspec(dllexport) int processBeamformingAndPostGPU(BeamformingGPUHandle handle, int8_t* input, int size, uint8_t* output, int* imagesize);

	//
	__declspec(dllexport) int processULMGPU(BeamformingGPUHandle handle, float* input_1, float* input_2, float* output, int step);

	//最后释放
	__declspec(dllexport) void deleteBeamformingGPUHandle(BeamformingGPUHandle handle);

	__declspec(dllexport) void deleteULMGPUHandle(BeamformingGPUHandle handle);
#ifdef __cplusplus
}
#endif
