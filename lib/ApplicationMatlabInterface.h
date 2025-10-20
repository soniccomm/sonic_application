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
		float element_radius,
		int sample_n,
		int bf_sample_num,
		int BeamN,
		float sos,
		float fs,
		float Js_Min_Aper,
		float Js_Max_Aper,
		float Js_Max_Fn,
		float Js_Min_Fn,
		float Js_Fn_depth,
		float Js_fn_setp,
		float Js_Fn_value,
		float* plae_wave_steer, int plae_wave_steer_size,
		float* beamx, int beamx_size,
		float* beamz, int beamz_size,
		float* cstartoffset, int cstartoffset_size,
		float* Jswin, int Jswin_size,
		float* channel_map, int channel_map_size,
		float* element_pos_xz, int element_pos_xz_size
	);

	__declspec(dllexport) BeamformingGPUHandle initializeBfiqGPU(
		int numperfile,
		int rx_num,
		int Channel,
		int RxChannel,
		int element_num,
		float element_pitch,
		float element_radius,
		int sample_n,
		int bf_sample_num,
		int BeamN,
		float sos,
		float fs,
		float fc,//发射频率
		float* plae_wave_steer, int plae_wave_steer_size,
		float* beamx, int beamx_size,
		float* beamz, int beamz_size,
		float* channel_map, int channel_map_size,
		float* element_pos_xz, int element_pos_xz_size,
		float* filter_coef, int filter_size,
		float* rx_apod, int rx_apod_size,
		float* rx_delay, int rx_delay_size,
		float* tx_delay, int tx_delay_size,
		float* angle_weight_table, int angle_weight_size
	);

	__declspec(dllexport) BeamformingGPUHandle initializespectralGPU(
		int prf,		int numperfile,
		int buffer_num,		int fft_period,		int lag,		
		float t_axis_span,		int t_idx,		int t_buffer_idx,		
		const char* probe_type,
		int Demod_AFE_Dynamic,
		int rx_num,
		int Channel,		int RxChannel,
		int element_num,		float element_pitch,
		int adc_head,
		int bf_sample_num,
		int BeamN,
		float sos,
		float fs,
		float fc, //发射频率
		float* plae_wave_steer, int plae_wave_steer_size,
		float* beamx, int beamx_size,
		float* beamz, int beamz_size,
		float* cstartoffset, int cstartoffset_size,
		float* channel_map, int channel_map_size,
		float* element_pos_xz, int element_pos_xz_size,
		float* filter_coef, int filter_size,
		float* rx_apod, int rx_apod_size,
		float* rx_delay, int rx_delay_size,
		float* tx_delay, int tx_delay_size,
		float* x_loc_all, float* z_loc_all, int loc_num,
		float* iir_b, float* iir_a, float* iir_zi, int iir_order
	);

	__declspec(dllexport) BeamformingGPUHandle initializecolorGPU(
		int prf, int numperfile,		int buffer_num, int fft_period, 
		int lag,		float t_axis_span, int t_idx, int t_buffer_idx,
		const char* probe_type,
		int Demod_AFE_Dynamic,
		int rx_num,
		int Channel, int RxChannel,
		int element_num, float element_pitch,
		int adc_head,
		int bf_sample_num,
		int BeamN,
		float sos,
		float fs,
		float fc, //发射频率
		int svd_auto,
		int svd_ord1,
		int svd_ord2,
		float* plae_wave_steer, int plae_wave_steer_size,
		float* beamx, int beamx_size,
		float* beamz, int beamz_size,
		float* cstartoffset, int cstartoffset_size,
		float* channel_map, int channel_map_size,
		float* element_pos_xz, int element_pos_xz_size,
		float* filter_coef, int filter_size,
		float* rx_apod, int rx_apod_size,
		float* rx_delay, int rx_delay_size,
		float* tx_delay, int tx_delay_size,
		float* x_loc_all, float* z_loc_all, int loc_num,
		int region_width, int region_height
	);

	__declspec(dllexport) BeamformingGPUHandle initializepowerGPU(
		int prf, int numperfile, int buffer_num, int fft_period,
		int lag, float t_axis_span, int t_idx, int t_buffer_idx,
		const char* probe_type,
		int Demod_AFE_Dynamic,
		int rx_num,
		int Channel, int RxChannel,
		int element_num, float element_pitch,
		int adc_head,
		int bf_sample_num,
		int BeamN,
		float sos,
		float fs,
		float fc, //发射频率
		int svd_auto,
		int svd_ord1,
		int svd_ord2,
		float* plae_wave_steer, int plae_wave_steer_size,
		float* beamx, int beamx_size,
		float* beamz, int beamz_size,
		float* cstartoffset, int cstartoffset_size,
		float* channel_map, int channel_map_size,
		float* element_pos_xz, int element_pos_xz_size,
		float* filter_coef, int filter_size,
		float* rx_apod, int rx_apod_size,
		float* rx_delay, int rx_delay_size,
		float* tx_delay, int tx_delay_size,
		float* x_loc_all, float* z_loc_all, int loc_num,
		int region_width, int region_height
	);

	 
	__declspec(dllexport) 	BeamformingGPUHandle  initializeAllGPU(
		const char* probe_type,	int Demod_AFE_Dynamic,	int rx_num,
		int Channel,int RxChannel,int element_num,	float element_pitch, float element_radius,
		int adc_head, int sample_n, int bf_sample_num,	int BeamN,
		float sos,	float fs,	float focus_depth,
		float Js_Min_Aper,	float Js_Max_Aper,	float Js_Max_Fn,	float Js_Min_Fn,	float Js_Fn_depth,	float Js_fn_setp,	float Js_Fn_value,
		float Fsweep_step,	float Fsweep_start,	float Fsweep_end,	float Roi_start,
		float* plae_wave_steer, int plae_wave_steer_size,
		float* beamx, int beamx_size,	float* beamz, int beamz_size,
		float* cstartoffset, int cstartoffset_size,
		float* Jswin, int Jswin_size,
		float* udm_L, int udm_L_size,
		float* fs_focus, int fs_focus_size,
		float* channel_map, int channel_map_size,
		float* element_pos_xz, int element_pos_xz_size,	
		int ValidBeamN,	int rx_line_size,	float* rx_line,	float SmoothCoef,	float UI_FrameRate, 	float B_fs,
		int Demod_MaxPoint,	float ImageStart,	float UI_depth,
		int dr_input,	int dr_output,	int dr_drange,
		int axialGain_size,	float* axialGain,	float xuanniuGain,	int spa_smoothcoef_size,	float* spa_smoothcoef,
		float* Demodulate_FreDepth,	float* Demodulate_FreValue,
		const char* Demodulate_FilterType,
		float blackholethre,
		float decifs,
		float* Demod_FilterDepth,	float* Demod_FilterCutoff,
		int  Demod_OrderFactor,
		float wininfo_alpha,
		int UI_DscHeight, int UI_DscWidth
	);
	
	__declspec(dllexport) 	BeamformingGPUHandle  initializeAllGPU_Test();
	//实时调用
	//return: <0  失败 >0 表示没有形成一张图,比如平面波的前面几帧  ==0 成功
	__declspec(dllexport) int processDataBeamformingGPU(BeamformingGPUHandle handle, int8_t* input, int size, float* output);

	__declspec(dllexport) int processDataBeamformingIQGPU(BeamformingGPUHandle handle, int8_t* input, int size, float* output);

	__declspec(dllexport) int processSpectralDataBeamformingandPostGPU(BeamformingGPUHandle handle, int8_t* input, 
		int size, float* output_b, float* output_sd, int bag_idx);

	__declspec(dllexport) int processColorDataBeamformingandPostGPU(BeamformingGPUHandle handle, int8_t* input,
		int size, float* output_b, float* output_cd, float* output_fm, int bag_idx);

	__declspec(dllexport) int processPowerDataBeamformingandPostGPU(BeamformingGPUHandle handle, int8_t* input,
		int size, float* output_b, float* output_cd, int bag_idx);
	
	//后处理
	__declspec(dllexport) int processBeamformingAndPostGPU(BeamformingGPUHandle handle, int8_t* input, int size, uint8_t* output, int* imagesize);

	//最后释放
	__declspec(dllexport) void deleteBeamformingGPUHandle(BeamformingGPUHandle handle);
#ifdef __cplusplus
}
#endif
