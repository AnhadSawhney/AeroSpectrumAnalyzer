#pragma once
struct Measure {
public: 
	enum Port {
		PORT_OUTPUT,
		PORT_INPUT,
	};
	enum Channel {
		CHANNEL_FL,
		CHANNEL_FR,
		CHANNEL_C,
		CHANNEL_LFE,
		CHANNEL_BL,
		CHANNEL_BR,
		CHANNEL_SL,
		CHANNEL_SR,
		MAX_CHANNELS,
		CHANNEL_SUM = MAX_CHANNELS
	};
	enum Type {
		TYPE_RMS,
		TYPE_PEAK,
		TYPE_FFT,
		TYPE_BAND,
		TYPE_FFTFREQ,
		TYPE_BANDFREQ,
		TYPE_FORMAT,
		TYPE_DEV_STATUS,
		TYPE_DEV_NAME,
		TYPE_DEV_ID,
		TYPE_DEV_LIST,
		// ... //
		NUM_TYPES
	};
	enum Format {
		FMT_INVALID,
		FMT_PCM_S16,
		FMT_PCM_F32,
		// ... //
		NUM_FORMATS
	};

	struct BandInfo {
		float	freq;
		float	x;
	};

	Port					m_port;						// port specifier (parsed from options)
	Channel					m_channel;					// channel specifier (parsed from options)
	//Type					m_type;						// data type specifier (parsed from options)
	Format					m_format;					// format specifier (detected in init)
	//int						m_envRMS[2];				// RMS attack/decay times in ms (parsed from options)
	//int						m_envPeak[2];				// peak attack/decay times in ms (parsed from options)
	int						m_envFFT[2];				// FFT attack/decay times in ms (parsed from options)
	int						m_fftSize;					// size of FFT (parsed from options)
	int						m_fftOverlap;				// number of samples between FFT calculations
	int						m_fftIdx;					// FFT index to retrieve (parsed from options)
	//int						m_nBands;					// number of frequency bands (parsed from options)
	//int						m_bandIdx;					// band index to retrieve (parsed from options)
	//double					m_gainRMS;					// RMS gain (parsed from options)
	//double					m_gainPeak;					// peak gain (parsed from options)
	double					m_freqMin;					// min freq for band measurement
	double					m_freqMax;					// max freq for band measurement
	double					m_sensitivity;				// dB range for FFT/Band return values (parsed from options)
	//Measure* m_parent;					// parent measure, if any
	//void* m_skin;						// skin pointer
	//LPCWSTR					m_rmName;					// measure name
	IMMDeviceEnumerator* m_enum;						// audio endpoint enumerator
	IMMDevice* m_dev;						// audio endpoint device
	WAVEFORMATEX* m_wfx;						// audio format info
	IAudioClient* m_clAudio;					// audio client instance
	IAudioCaptureClient* m_clCapture;				// capture client instance
#if (WINDOWS_BUG_WORKAROUND)
	IAudioClient* m_clBugAudio;				// audio client for dummy silent channel
	IAudioRenderClient* m_clBugRender;				// render client for dummy silent channel
#endif
	WCHAR					m_reqID[64];				// requested device ID (parsed from options)
	WCHAR					m_devName[64];				// device friendly name (detected in init)
	//float					m_kRMS[2];					// RMS attack/decay filter constants
	//float					m_kPeak[2];					// peak attack/decay filter constants
	float					m_kFFT[2];					// FFT attack/decay filter constants
	//double					m_rms[MAX_CHANNELS];		// current RMS levels
	//double					m_peak[MAX_CHANNELS];		// current peak levels
	double					m_pcMult;					// performance counter inv frequency
	LARGE_INTEGER			m_pcFill;					// performance counter on last full buffer
	LARGE_INTEGER			m_pcPoll;					// performance counter on last device poll
	kiss_fftr_cfg			m_fftCfg;// [MAX_CHANNELS] ;		// FFT states for each channel
	float* m_fftIn;// [MAX_CHANNELS] ;		// buffer for each channel's FFT input
	float* m_fftOut;// [MAX_CHANNELS] ;		// buffer for each channel's FFT output
	float* m_fftKWdw;					// window function coefficients
	float* m_fftTmpIn;					// temp FFT processing buffer
	kiss_fft_cpx* m_fftTmpOut;				// temp FFT processing buffer
	int						m_fftBufW;					// write index for input ring buffers
	int						m_fftBufP;					// decremental counter - process FFT at zero
	//float* m_bandFreq;					// buffer of band max frequencies
	//float* m_bandOut; // [MAX_CHANNELS];	// buffer of band values
	float m_bandFreq[BANDS];
	float m_bandOut[BANDS];

	Measure(Options* o) :
		m_port(PORT_OUTPUT),
		m_channel(CHANNEL_SUM),
		//m_type(TYPE_RMS),
		m_format(FMT_INVALID),
		m_fftSize(0),
		m_fftOverlap(0),
		m_fftIdx(-1),
		//m_nBands(0),
		//m_bandIdx(-1),
		//m_gainRMS(1.0),
		//m_gainPeak(1.0),
		m_freqMin(20.0),
		m_freqMax(20000.0),
		m_sensitivity(35.0),
		//m_parent(NULL),
		//m_skin(NULL),
		//m_rmName(NULL),
		m_enum(NULL),
		m_dev(NULL),
		m_wfx(NULL),
		m_clAudio(NULL),
		m_clCapture(NULL),
#if (WINDOWS_BUG_WORKAROUND)
		m_clBugAudio(NULL),
		m_clBugRender(NULL),
#endif
		m_fftKWdw(NULL),
		m_fftTmpIn(NULL),
		m_fftTmpOut(NULL),
		m_fftBufW(0),
		m_fftBufP(0)
		//m_bandFreq(NULL)
	{
		//m_envRMS[0] = 300;
		//m_envRMS[1] = 300;
		//m_envPeak[0] = 50;
		//m_envPeak[1] = 2500;
		m_envFFT[0] = 300;
		m_envFFT[1] = 300;
		m_reqID[0] = '\0';
		m_devName[0] = '\0';
		//m_kRMS[0] = 0.0f;
		//m_kRMS[1] = 0.0f;
		//m_kPeak[0] = 0.0f;
		//m_kPeak[1] = 0.0f;
		m_kFFT[0] = 0.0f;
		m_kFFT[1] = 0.0f;
		//for (int iChan = 0; iChan < MAX_CHANNELS; ++iChan) {
			//m_rms[iChan] = 0.0;
			//m_peak[iChan] = 0.0;
			m_fftCfg/*[iChan]*/ = NULL;
			m_fftIn/*[iChan]*/ = NULL;
			m_fftOut/*[iChan]*/ = NULL;
		for (int i = 0; i < BANDS; i++) {
			m_bandOut[i] = NULL;
			m_bandFreq[i] = NULL;
		}
		//}
		LARGE_INTEGER pcFreq;
		QueryPerformanceFrequency(&pcFreq);
		m_pcMult = 1.0 / (double)pcFreq.QuadPart;
		Initialize(o);
	}

	void Initialize(Options* o);
	void Finalize();
	void Reload(Options* o);
	void Update();

	HRESULT	DeviceInit();
	void DeviceRelease();
};

const CLSID CLSID_MMDeviceEnumerator = __uuidof(MMDeviceEnumerator);
const IID IID_IMMDeviceEnumerator = __uuidof(IMMDeviceEnumerator);
const IID IID_IAudioClient = __uuidof(IAudioClient);
const IID IID_IAudioCaptureClient = __uuidof(IAudioCaptureClient);
const IID IID_IAudioRenderClient = __uuidof(IAudioRenderClient);

/**
 * Create and initialize a measure instance.  Creates WASAPI loopback
 * device if not a child measure.
 *
 * @param[out]	data			Pointer address in which to return measure instance.
 * @param[in]	rm				Rainmeter context.
 */
void Measure::Initialize(Options* o) {
	//m_skin = RmGetSkin(rm);
	//m_rmName = RmGetMeasureName(rm);

	// parse port specifier
	LPCWSTR port = o->readLPCWSTR("Port");
	if (port && *port) {
		if (_wcsicmp(port, L"Output") == 0) {
			m_port = Measure::PORT_OUTPUT;
		} else if (_wcsicmp(port, L"Input") == 0) {
			m_port = Measure::PORT_INPUT;
		} else {
			std::cerr << "ERROR: Invalid Port " << port <<" must be one of: Output or Input.\n";
		}
	}

	// parse requested device ID (optional)
	LPCWSTR reqID = o->readLPCWSTR("ID");
	if (reqID) {
		_snwprintf_s(m_reqID, _TRUNCATE, L"%s", reqID);
	}

	// initialize FFT data
	m_fftSize = o->readInt("Size");
	if (m_fftSize < 0 || (m_fftSize & 1)) {
		std::cerr << "ERROR: Invalid FFTSize " << m_fftSize << " must be an even integer >= 0. (powers of 2 work best)\n";
		m_fftSize = 0;
	}

	if (m_fftSize) {
		m_fftOverlap = o->readInt("Overlap");
		if(m_fftOverlap < 0	|| m_fftOverlap >= m_fftSize) {
			std::cerr << "ERROR: Invalid FFTOverlap " << m_fftOverlap << ": must be an integer between 0 and FFTSize (" << m_fftSize << ").\n";
			m_fftOverlap = 0;
		}
	}

	// initialize frequency bands
	/*m_nBands = o->readInt("Bands");
	if (m_nBands < 0) {
		std::cerr << "ERROR: Invalid Bands " << m_nBands << ": must be an integer >= 0.\n";
		m_nBands = 0;
	}*/
	m_freqMin = max(0.0, o->readDouble("Minimum Frequency"));
	m_freqMax = max(0.0, o->readDouble("Maximum Frequency"));

	// initialize the watchdog timer
	QueryPerformanceCounter(&m_pcPoll);

	// create the enumerator
	if (CoCreateInstance(CLSID_MMDeviceEnumerator, NULL, CLSCTX_ALL, IID_IMMDeviceEnumerator, (void**)&m_enum) == S_OK) {
		// init the device (ok if it fails - it'll keep checking during Update)
		DeviceInit();
		return;
	}

	Reload(o);

	SAFE_RELEASE(m_enum)
}


/**
 * Destroy the measure instance.
 *
 * @param[in]	data			Measure instance pointer.
 */
void Measure::Finalize() {
	DeviceRelease();
	SAFE_RELEASE(m_enum)
	delete this;
}

/**
 * (Re-)parse parameters from .ini file.
 *
 * @param[in]	data			Measure instance pointer.
 * @param[in]	rm				Rainmeter context.
 * @param[out]	maxValue		?
 */
void Measure::Reload(Options* o) {
	/*static const LPCWSTR s_typeName[Measure::NUM_TYPES] = {
		L"RMS",								// TYPE_RMS
		L"Peak",							// TYPE_PEAK
		L"FFT",								// TYPE_FFT
		L"Band",							// TYPE_BAND
		L"FFTFreq",							// TYPE_FFTFREQ
		L"BandFreq",						// TYPE_BANDFREQ
		L"Format",							// TYPE_FORMAT
		L"DeviceStatus",					// TYPE_DEV_STATUS
		L"DeviceName",						// TYPE_DEV_NAME
		L"DeviceID",						// TYPE_DEV_ID
		L"DeviceList",						// TYPE_DEV_LIST
	};*/
	static const LPCWSTR s_chanName[Measure::CHANNEL_SUM + 1][3] = {
		{ L"L",		L"FL",		L"0", },	// CHANNEL_FL
		{ L"R",		L"FR",		L"1", },	// CHANNEL_FR
		{ L"C",		L"",		L"2", },	// CHANNEL_C
		{ L"LFE",	L"Sub",		L"3", },	// CHANNEL_LFE
		{ L"BL",	L"",		L"4", },	// CHANNEL_BL
		{ L"BR",	L"",		L"5", },	// CHANNEL_BR
		{ L"SL",	L"",		L"6", },	// CHANNEL_SL
		{ L"SR",	L"",		L"7", },	// CHANNEL_SR
		{ L"Sum",	L"Avg",		L"", },		// CHANNEL_SUM
	};

	// parse channel specifier
	LPCWSTR channel = o->readLPCWSTR("Channel");
	if (*channel) {
		bool found = false;
		for (int iChan = 0; iChan <= Measure::CHANNEL_SUM && !found; ++iChan) {
			for (int j = 0; j < 3; ++j) {
				if (_wcsicmp(channel, s_chanName[iChan][j]) == 0) {
					m_channel = (Measure::Channel)iChan;
					found = true;
					break;
				}
			}
		}
		if (!found) {
			WCHAR	msg[512];
			WCHAR* d = msg;
			d += _snwprintf_s(d, (sizeof(msg) + (UINT32)msg - (UINT32)d) / sizeof(WCHAR), _TRUNCATE, L"Invalid Channel '%s', must be an integer between 0 and %d, or one of:", channel,	Measure::MAX_CHANNELS - 1);
			for (unsigned int i = 0; i <= Measure::CHANNEL_SUM; ++i) {
				d += _snwprintf_s(d, (sizeof(msg) + (UINT32)msg - (UINT32)d) / sizeof(WCHAR), _TRUNCATE, (i? L", %s%s" : L" %s%s"),	(i == Measure::CHANNEL_SUM) ? L"or " : L"",	s_chanName[i][0]);
			}
			d += _snwprintf_s(d, (sizeof(msg) + (UINT32)msg - (UINT32)d) / sizeof(WCHAR), _TRUNCATE, L".\n");
			std::cerr << "ERROR: " << msg;
		}
	}
// parse data type
/*LPCWSTR type = RmReadString(rm, L"Type", L"");
if (*type) {
	int iType;
	for (iType = 0; iType < Measure::NUM_TYPES; ++iType) {
		if (_wcsicmp(type, s_typeName[iType]) == 0) {
			m_type = (Measure::Type)iType;
			break;
		}
	}
	if (!(iType < Measure::NUM_TYPES)) {
		WCHAR	msg[512];
		WCHAR* d = msg;
		d += _snwprintf_s(d, (sizeof(msg) + (UINT32)msg - (UINT32)d) / sizeof(WCHAR), _TRUNCATE, L"Invalid Type '%s', must be one of:", type);
		for (unsigned int i = 0; i < Measure::NUM_TYPES; ++i) {
			d += _snwprintf_s(d, (sizeof(msg) + (UINT32)msg - (UINT32)d) / sizeof(WCHAR), _TRUNCATE, (i ? L", %s%s" : L" %s%s"), (i == (Measure::NUM_TYPES - 1)) ? L"or " : L"", s_typeName[i]);
		}
		d += _snwprintf_s(d, (sizeof(msg) + (UINT32)msg - (UINT32)d) / sizeof(WCHAR), _TRUNCATE, L".\n");
		ERROR_LOG("AudioLevel", msg);
	}
}
// parse FFT index request
m_fftIdx = max(0, RmReadInt(rm, L"FFTIdx", m_fftIdx));
if (m_parent) {
	m_fftIdx = min(m_parent->m_fftSize / 2, m_fftIdx);
}
else {
	m_fftIdx = min(m_fftSize / 2, m_fftIdx);
}
// parse band index request
m_bandIdx = max(0, RmReadInt(rm, L"BandIdx", m_bandIdx));
if (m_parent) {
	m_bandIdx = min(m_parent->m_nBands, m_bandIdx);
}
else {
	m_bandIdx = min(m_nBands, m_bandIdx);
}*/

// parse envelope values on parents only
//if (!m_parent) {
	// (re)parse envelope values
//	m_envRMS[0] = max(0, RmReadInt(rm, L"RMSAttack", m_envRMS[0]));
//	m_envRMS[1] = max(0, RmReadInt(rm, L"RMSDecay", m_envRMS[1]));
//	m_envPeak[0] = max(0, RmReadInt(rm, L"PeakAttack", m_envPeak[0]));
//	m_envPeak[1] = max(0, RmReadInt(rm, L"PeakDecay", m_envPeak[1]));
	this -> m_envFFT[0] = max(0, o->readInt("Attack"));
	this -> m_envFFT[1] = max(0, o->readInt("Decay"));
	// (re)parse gain constants
//	m_gainRMS = max(0.0, RmReadDouble(rm, L"RMSGain", m_gainRMS));
//	m_gainPeak = max(0.0, RmReadDouble(rm, L"PeakGain", m_gainPeak));
	m_sensitivity = max(1.0, o->readDouble("Sensitivity"));
	// regenerate filter constants
	if (m_wfx) {
		const double freq = m_wfx->nSamplesPerSec;
//		m_kRMS[0] = (float)exp(log10(0.01) / (freq * (double)m_envRMS[0] * 0.001));
//		m_kRMS[1] = (float)exp(log10(0.01) / (freq * (double)m_envRMS[1] * 0.001));
//		m_kPeak[0] = (float)exp(log10(0.01) / (freq * (double)m_envPeak[0] * 0.001));
//		m_kPeak[1] = (float)exp(log10(0.01) / (freq * (double)m_envPeak[1] * 0.001));
//		if (m_fftSize) {
			m_kFFT[0] = (float)exp(log10(0.01) / (freq / (m_fftSize - m_fftOverlap) * (double)m_envFFT[0] * 0.001));
			m_kFFT[1] = (float)exp(log10(0.01) / (freq / (m_fftSize - m_fftOverlap) * (double)m_envFFT[1] * 0.001));
//		}
//	}
	}
}


/**
 * Update the measure.
 *
 * @param[in]	data			Measure instance pointer.
 * @return		Latest value - typically an audio level between 0.0 and 1.0.
 */
void Measure::Update() {
	LARGE_INTEGER	pcCur; 
	QueryPerformanceCounter(&pcCur);

	// query the buffer
	if (m_clCapture	&& ((pcCur.QuadPart - m_pcPoll.QuadPart) * m_pcMult) >= QUERY_TIMEOUT) {
		BYTE* buffer;
		UINT32			nFrames;
		DWORD			flags;
		UINT64			pos;
		HRESULT			hr;
		while ((hr = m_clCapture->GetBuffer(&buffer, &nFrames, &flags, &pos, NULL)) == S_OK) {
			// measure RMS and peak levels
			/*float		rms[Measure::MAX_CHANNELS];
			float		peak[Measure::MAX_CHANNELS];
			for (int iChan = 0; iChan < Measure::MAX_CHANNELS; ++iChan) {
				rms[iChan] = (float)m_rms[iChan];
				peak[iChan] = (float)m_peak[iChan];
			}
			// loops unrolled for float, 16b and mono, stereo
			if (m_format == Measure::FMT_PCM_F32) {
				float* s = (float*)buffer;
				if (m_wfx->nChannels == 1) {
					for (unsigned int iFrame = 0; iFrame < nFrames; ++iFrame) {
						float xL = (float)*s++;
						float sqrL = xL * xL;
						float absL = abs(xL);
						rms[0] = sqrL + m_kRMS[(sqrL < rms[0])] * (rms[0] - sqrL);
						peak[0] = absL + m_kPeak[(absL < peak[0])] * (peak[0] - absL);
						rms[1] = rms[0];
						peak[1] = peak[0];
					}
				} else if (m_wfx->nChannels == 2) {
					for (unsigned int iFrame = 0; iFrame < nFrames; ++iFrame) {
						float xL = (float)*s++;
						float xR = (float)*s++;
						float sqrL = xL * xL;
						float sqrR = xR * xR;
						float absL = abs(xL);
						float absR = abs(xR);
						rms[0] = sqrL + m_kRMS[(sqrL < rms[0])] * (rms[0] - sqrL);
						rms[1] = sqrR + m_kRMS[(sqrR < rms[1])] * (rms[1] - sqrR);
						peak[0] = absL + m_kPeak[(absL < peak[0])] * (peak[0] - absL);
						peak[1] = absR + m_kPeak[(absR < peak[1])] * (peak[1] - absR);
					}
				} else {
					for (unsigned int iFrame = 0; iFrame < nFrames; ++iFrame) {
						for (unsigned int iChan = 0; iChan < m_wfx->nChannels; ++iChan) {
							float	x = (float)*s++;
							float	sqrX = x * x;
							float	absX = abs(x);
							rms[iChan] = sqrX + m_kRMS[(sqrX < rms[iChan])] * (rms[iChan] - sqrX);
							peak[iChan] = absX + m_kPeak[(absX < peak[iChan])] * (peak[iChan] - absX);
						}
					}
				}
			} else if (m_format == Measure::FMT_PCM_S16) {
				INT16* s = (INT16*)buffer;
				if (m_wfx->nChannels == 1) {
					for (unsigned int iFrame = 0; iFrame < nFrames; ++iFrame) {
						float xL = (float)*s++ * 1.0f / 0x7fff;
						float sqrL = xL * xL;
						float absL = abs(xL);
						rms[0] = sqrL + m_kRMS[(sqrL < rms[0])] * (rms[0] - sqrL);
						peak[0] = absL + m_kPeak[(absL < peak[0])] * (peak[0] - absL);
						rms[1] = rms[0];
						peak[1] = peak[0];
					}
				} else if (m_wfx->nChannels == 2) {
					for (unsigned int iFrame = 0; iFrame < nFrames; ++iFrame) {
						float xL = (float)*s++ * 1.0f / 0x7fff;
						float xR = (float)*s++ * 1.0f / 0x7fff;
						float sqrL = xL * xL;
						float sqrR = xR * xR;
						float absL = abs(xL);
						float absR = abs(xR);
						rms[0] = sqrL + m_kRMS[(sqrL < rms[0])] * (rms[0] - sqrL);
						rms[1] = sqrR + m_kRMS[(sqrR < rms[1])] * (rms[1] - sqrR);
						peak[0] = absL + m_kPeak[(absL < peak[0])] * (peak[0] - absL);
						peak[1] = absR + m_kPeak[(absR < peak[1])] * (peak[1] - absR);
					}
				} else {
					for (unsigned int iFrame = 0; iFrame < nFrames; ++iFrame) {
						for (unsigned int iChan = 0; iChan < m_wfx->nChannels; ++iChan) {
							float	x = (float)*s++ * 1.0f / 0x7fff;
							float	sqrX = x * x;
							float	absX = abs(x);
							rms[iChan] = sqrX + m_kRMS[(sqrX < rms[iChan])] * (rms[iChan] - sqrX);
							peak[iChan] = absX + m_kPeak[(absX < peak[iChan])] * (peak[iChan] - absX);
						}
					}
				}
			}
			for (int iChan = 0; iChan < Measure::MAX_CHANNELS; ++iChan) {
				m_rms[iChan] = rms[iChan];
				m_peak[iChan] = peak[iChan];
			}*/

			// process FFTs (optional)
			float* sF32 = (float*)buffer;
			INT16* sI16 = (INT16*)buffer;
			const float	scalar = (float)(1.0 / sqrt(m_fftSize));
			for (unsigned int iFrame = 0; iFrame < nFrames; ++iFrame) {
				// fill ring buffers (demux streams)
				//for (unsigned int iChan = 0; iChan < m_wfx->nChannels; ++iChan) {
					(m_fftIn/*[iChan]*/)[m_fftBufW] = (m_format == Measure::FMT_PCM_F32) ? *sF32++ : ((float)*sI16++ * 1.0f / 0x7fff);
				//}
				m_fftBufW = (m_fftBufW + 1) % m_fftSize;
				// if overlap limit reached, process FFTs for each channel
				if (!--m_fftBufP) {
					//for (unsigned int iChan = 0; iChan < m_wfx->nChannels; ++iChan) {
						if (!(flags & AUDCLNT_BUFFERFLAGS_SILENT)) {
							// copy from the ring buffer to temp space
							memcpy(&m_fftTmpIn[0], &(m_fftIn/*[iChan]*/)[m_fftBufW], (m_fftSize - m_fftBufW) * sizeof(float));
							memcpy(&m_fftTmpIn[m_fftSize - m_fftBufW], &m_fftIn/*[iChan]*/[0], m_fftBufW * sizeof(float));
							// apply the windowing function
							for (int iBin = 0; iBin < m_fftSize; ++iBin) {
								m_fftTmpIn[iBin] *= m_fftKWdw[iBin];
							}
							kiss_fftr(m_fftCfg/*[iChan]*/, m_fftTmpIn, m_fftTmpOut);
						}
						else {
							memset(m_fftTmpOut, 0, m_fftSize * sizeof(kiss_fft_cpx));
						}
						// filter the bin levels as with peak measurements
						for (int iBin = 0; iBin < m_fftSize; ++iBin) {
							float x0 = (m_fftOut/*[iChan]*/)[iBin];
							float x1 = (m_fftTmpOut[iBin].r * m_fftTmpOut[iBin].r + m_fftTmpOut[iBin].i * m_fftTmpOut[iBin].i) * scalar;
							x0 = x1 + m_kFFT[(x1 < x0)] * (x0 - x1);
							(m_fftOut/*[iChan]*/)[iBin] = x0;
						}
					//}
					m_fftBufP = m_fftSize - m_fftOverlap;
				}

				// integrate FFT results into log-scale frequency bands
				const float		df = (float)m_wfx->nSamplesPerSec / m_fftSize;
				const float		scalar = 2.0f / (float)m_wfx->nSamplesPerSec;
				//for (unsigned int iChan = 0; iChan < m_wfx->nChannels; ++iChan) {
					memset(m_bandOut/*[iChan]*/, 0, BANDS * sizeof(float));
					int			iBin = 0;
					int			iBand = 0;
					float		f0 = 0.0f;
					while (iBin <= (m_fftSize / 2) && iBand < BANDS) {
						float	fLin1 = ((float)iBin + 0.5f) * df;
						float	fLog1 = m_bandFreq[iBand];
						float	x = (m_fftOut/*[iChan]*/)[iBin];
						float& y = (m_bandOut/*[iChan]*/)[iBand];
						if (fLin1 <= fLog1) {
							y += (fLin1 - f0) * x * scalar;
							f0 = fLin1;
							iBin += 1;
						}
						else {
							y += (fLog1 - f0) * x * scalar;
							f0 = fLog1;
							iBand += 1;
						}
					}
				//}

				// release the buffer
				m_clCapture->ReleaseBuffer(nFrames);

				// mark the time of last buffer update
				m_pcFill = pcCur;
			}
		}
		// detect device disconnection
		switch (hr) {
			case AUDCLNT_S_BUFFER_EMPTY:
				// Windows bug: sometimes when shutting down a playback application, it doesn't zero
				// out the buffer.  Detect this by checking the time since the last successful fill
				// and resetting the volumes if past the threshold.
				/*if (((pcCur.QuadPart - m_pcFill.QuadPart) * m_pcMult) >= EMPTY_TIMEOUT) {
					for (int iChan = 0; iChan < Measure::MAX_CHANNELS; ++iChan) {
						m_rms[iChan] = 0.0;
						m_peak[iChan] = 0.0;
					}
				}*/
				break;
			case AUDCLNT_E_BUFFER_ERROR:
			case AUDCLNT_E_DEVICE_INVALIDATED:
			case AUDCLNT_E_SERVICE_NOT_RUNNING:
				DeviceRelease();
				break;
		}
		m_pcPoll = pcCur;
	} else if (!m_clCapture	&& (((pcCur.QuadPart - m_pcPoll.QuadPart) * m_pcMult) >= DEVICE_TIMEOUT)) {
		// poll for new devices
		assert(m_enum);
		assert(!m_dev);
		DeviceInit();
		m_pcPoll = pcCur;
	}

	/*switch (m_type) {
		case Measure::TYPE_FFT:
			if (parent->m_clCapture && parent->m_fftSize) {
				double	x;
				const int iFFT = m_fftIdx;
				if (m_channel == Measure::CHANNEL_SUM) {
					if (parent->m_wfx->nChannels >= 2) {
						x = (parent->m_fftOut[0][iFFT] + parent->m_fftOut[1][iFFT]) * 0.5;
					} else {
						x = parent->m_fftOut[0][iFFT];
					}
				} else if (m_channel < parent->m_wfx->nChannels) {
					x = parent->m_fftOut[m_channel][iFFT];
				}
				x = CLAMP01(x);
				x = max(0, 10.0 / parent->m_sensitivity * log10(x) + 1.0);
				return x;
			}
			break;
		case Measure::TYPE_BAND:
			if (parent->m_clCapture && parent->m_nBands) {
				double	x;
				const int iBand = m_bandIdx;
				if (m_channel == Measure::CHANNEL_SUM) {
					if (parent->m_wfx->nChannels >= 2) {
						x = (parent->m_bandOut[0][iBand] + parent->m_bandOut[1][iBand]) * 0.5;
					} else {
						x = parent->m_bandOut[0][iBand];
					}
				} else if (m_channel < parent->m_wfx->nChannels) {
					x = parent->m_bandOut[m_channel][iBand];
				}
				x = CLAMP01(x);
				x = max(0, 10.0 / parent->m_sensitivity * log10(x) + 1.0);
				return x;
			}
			break;
		case Measure::TYPE_FFTFREQ:
			if (parent->m_clCapture && parent->m_fftSize && m_fftIdx <= (parent->m_fftSize / 2)) {
				return (m_fftIdx * m_wfx->nSamplesPerSec / parent->m_fftSize);
			}
			break;
		case Measure::TYPE_BANDFREQ:
			if (parent->m_clCapture	&& parent->m_nBands	&& m_bandIdx < (parent->m_nBands)) {
				return parent->m_bandFreq[m_bandIdx];
			}
			break;
		case Measure::TYPE_DEV_STATUS:
			if (parent->m_dev) {
				DWORD state;
				if (parent->m_dev->GetState(&state) == S_OK	&& state == DEVICE_STATE_ACTIVE) {
					return 1.0;
				}
			}
			break;
	}*/
	//return 0.0;
}

/**
 * Get a string value from the measure.
 *
 * @param[in]	data			Measure instance pointer.
 * @return		String value - must be copied out by the caller.
 */
/*LPCWSTR Measure::GetString() {
	//Measure* parent = (m_parent) ? m_parent : m;
	static WCHAR	buffer[512];
	const char* s_fmtName[Measure::NUM_FORMATS] = {
		"<invalid>",	// FMT_INVALID
		"PCM 16b",		// FMT_PCM_S16
		"PCM 32b",		// FMT_PCM_F32
	};

	buffer[0] = '\0';
	switch (m_type) {
		default:
			// return NULL for any numeric values, so Rainmeter can auto-convert them.
			return NULL;
		case Measure::TYPE_FORMAT:
			if (parent->m_wfx) {
				_snwprintf_s(buffer, _TRUNCATE,	L"%dHz %s %dch", parent->m_wfx->nSamplesPerSec,	s_fmtName[parent->m_format], parent->m_wfx->nChannels);
			}
			break;
		case Measure::TYPE_DEV_NAME:
			return parent->m_devName;
		case Measure::TYPE_DEV_ID:
			if (parent->m_dev) {
				LPWSTR	pwszID = NULL;
				if (parent->m_dev->GetId(&pwszID) == S_OK) {
					_snwprintf_s(buffer, _TRUNCATE, L"%s", pwszID);
					CoTaskMemFree(pwszID);
				}
			}
			break;
		case Measure::TYPE_DEV_LIST:
			if (parent->m_enum) {
				IMMDeviceCollection* collection = NULL;
				if (parent->m_enuEnumAudioEndpoints((parent->m_port == Measure::PORT_OUTPUT) ? eRender : eCapture, DEVICE_STATE_ACTIVE | DEVICE_STATE_UNPLUGGED, &collection) == S_OK) {
					WCHAR* d = &buffer[0];
					UINT nDevices;	
					collection->GetCount(&nDevices);
					for (ULONG iDevice = 0; iDevice < nDevices; ++iDevice) {
						IMMDevice* device = NULL;
						IPropertyStore* props = NULL;
						if (collection->Item(iDevice, &device) == S_OK && device->OpenPropertyStore(STGM_READ, &props) == S_OK) {
							LPWSTR		id = NULL;
							PROPVARIANT	varName;	PropVariantInit(&varName);
							if (device->GetId(&id) == S_OK && props->GetValue(PKEY_Device_FriendlyName, &varName) == S_OK) {
								d += _snwprintf_s(d, (sizeof(buffer) + (UINT32)buffer - (UINT32)d) / sizeof(WCHAR), _TRUNCATE, L"%s%s: %s", (iDevice > 0) ? "\n" : "", id, varName.pwszVal);
							}
							if (id) {
								CoTaskMemFree(id);
							}
							PropVariantClear(&varName);
						}
						SAFE_RELEASE(props);
						SAFE_RELEASE(device)
					}
				}
				SAFE_RELEASE(collection);
			}
			break;
	}
	return buffer;
}*/


/**
 * Try to initialize the default device for the specified port.
 *
 * @return		Result value, S_OK on success.
 */
HRESULT	Measure::DeviceInit() {
	HRESULT hr;

	// get the device handle
	assert(m_enum && !m_dev);

	// if a specific ID was requested, search for that one, otherwise get the default
	if (*m_reqID) {
		hr = m_enum->GetDevice(m_reqID, &m_dev);
		if (hr != S_OK)	{
			WCHAR msg[256];
			_snwprintf_s(msg, _TRUNCATE, L"Audio %s device '%s' not found (error 0x%08x).",	m_port == PORT_OUTPUT ? L"output" : L"input", m_reqID, hr);

			std::cerr << "WARNING: " << msg;
		}
	} else {
		hr = m_enum->GetDefaultAudioEndpoint(m_port == PORT_OUTPUT ? eRender : eCapture, eConsole, &m_dev);
	}

	

	// store device name
	IPropertyStore* props = NULL;

	double step;
	
	EXIT_ON_ERROR(hr);
	if (m_dev->OpenPropertyStore(STGM_READ, &props) == S_OK) {
		PROPVARIANT	varName;
		PropVariantInit(&varName);

		if (props->GetValue(PKEY_Device_FriendlyName, &varName) == S_OK) {
			_snwprintf_s(m_devName, _TRUNCATE, L"%s", varName.pwszVal);
		}

		PropVariantClear(&varName);
	}

	SAFE_RELEASE(props);

#if (WINDOWS_BUG_WORKAROUND)
	// get an extra audio client for the dummy silent channel
	hr = m_dev->Activate(IID_IAudioClient, CLSCTX_ALL, NULL, (void**)&m_clBugAudio);
	if (hr != S_OK)	{
		std::cerr << "WARNING: Failed to create audio client for Windows bug workaround.\n";
	}
#endif

	// get the main audio client
	hr = m_dev->Activate(IID_IAudioClient, CLSCTX_ALL, NULL, (void**)&m_clAudio);
	if (hr != S_OK)	{
		std::cerr << "WARNING: Failed to create audio client.\n";
	}

	EXIT_ON_ERROR(hr);

	// parse audio format - Note: not all formats are supported.
	hr = m_clAudio->GetMixFormat(&m_wfx);
	EXIT_ON_ERROR(hr);

	switch (m_wfx->wFormatTag)	{
		case WAVE_FORMAT_PCM:
			if (m_wfx->wBitsPerSample == 16) {
				m_format = FMT_PCM_S16;
			}
			break;

		case WAVE_FORMAT_IEEE_FLOAT:
			m_format = FMT_PCM_F32;
			break;

		case WAVE_FORMAT_EXTENSIBLE:
			if (reinterpret_cast<WAVEFORMATEXTENSIBLE*>(m_wfx)->SubFormat == KSDATAFORMAT_SUBTYPE_IEEE_FLOAT) {
				m_format = FMT_PCM_F32;
			}
			break;
	}

	if (m_format == FMT_INVALID) {
		std::cerr << "WARNING: Invalid sample format.  Only PCM 16b integer or PCM 32b float are supported.\n";
	}

	// setup FFT buffers
	//for (int iChan = 0; iChan < m_wfx->nChannels; ++iChan) {
		m_fftCfg/*[iChan]*/ = kiss_fftr_alloc(m_fftSize, 0, NULL, NULL);
		m_fftIn/*[iChan]*/ = (float*)calloc(m_fftSize * sizeof(float), 1);
		m_fftOut/*[iChan]*/ = (float*)calloc(m_fftSize * sizeof(float), 1);
	//}
	m_fftKWdw = (float*)calloc(m_fftSize * sizeof(float), 1);
	m_fftTmpIn = (float*)calloc(m_fftSize * sizeof(float), 1);
	m_fftTmpOut = (kiss_fft_cpx*)calloc(m_fftSize * sizeof(kiss_fft_cpx), 1);
	m_fftBufP = m_fftSize - m_fftOverlap;
	// calculate window function coefficients (http://en.wikipedia.org/wiki/Window_function#Hann_.28Hanning.29_window)
	for (int iBin = 0; iBin < m_fftSize; ++iBin) {
		m_fftKWdw[iBin] = (float)(0.5 * (1.0 - cos(TWOPI * iBin / (m_fftSize - 1))));
	}
	
	// calculate band frequencies and allocate band output buffers
	//m_bandFreq = (float*)malloc(BANDS * sizeof(float));
	step = (log(m_freqMax / m_freqMin) / BANDS) / log(2.0);
	m_bandFreq[0] = (float)(m_freqMin * pow(2.0, step / 2.0));

	for (int iBand = 1; iBand < BANDS; ++iBand) {
		m_bandFreq[iBand] = (float)(m_bandFreq[iBand - 1] * pow(2.0, step));
	}

	//for (int iChan = 0; iChan < m_wfx->nChannels; ++iChan) {
		//m_bandOut/*[iChan]*/ = (float*)calloc(BANDS * sizeof(float), 1);
	//}

	//REFERENCE_TIME hnsRequestedDuration = REFTIMES_PER_SEC;

#if (WINDOWS_BUG_WORKAROUND)
	// ---------------------------------------------------------------------------------------
	// Windows bug workaround: create a silent render client before initializing loopback mode
	// see: http://social.msdn.microsoft.com/Forums/windowsdesktop/en-US/c7ba0a04-46ce-43ff-ad15-ce8932c00171/loopback-recording-causes-digital-stuttering?forum=windowspro-audiodevelopment
	if (m_port == PORT_OUTPUT) {
		hr = m_clBugAudio->Initialize(AUDCLNT_SHAREMODE_SHARED, 0, REFTIMES_PER_SEC/*hnsRequestedDuration*/, 0, m_wfx, NULL);
		EXIT_ON_ERROR(hr);

		// get the frame count
		UINT32 nFrames;
		hr = m_clBugAudio->GetBufferSize(&nFrames);
		EXIT_ON_ERROR(hr);

		// create a render client
		hr = m_clBugAudio->GetService(IID_IAudioRenderClient, (void**)&m_clBugRender);
		EXIT_ON_ERROR(hr);

		// get the buffer
		BYTE* buffer;
		hr = m_clBugRender->GetBuffer(nFrames, &buffer);
		EXIT_ON_ERROR(hr);

		// release it
		hr = m_clBugRender->ReleaseBuffer(nFrames, AUDCLNT_BUFFERFLAGS_SILENT);
		EXIT_ON_ERROR(hr);

		// start the stream
		hr = m_clBugAudio->Start();
		EXIT_ON_ERROR(hr);
	}
	// ---------------------------------------------------------------------------------------
#endif

	// initialize the audio client
	hr = m_clAudio->Initialize(AUDCLNT_SHAREMODE_SHARED, m_port == PORT_OUTPUT ? AUDCLNT_STREAMFLAGS_LOOPBACK : 0, REFTIMES_PER_SEC/*hnsRequestedDuration*/, 0, m_wfx, NULL);
	if (hr != S_OK)	{
		// Compatibility with the Nahimic audio driver
		// https://github.com/rainmeter/rainmeter/commit/0a3dfa35357270512ec4a3c722674b67bff541d6
		// https://social.msdn.microsoft.com/Forums/windowsdesktop/en-US/bd8cd9f2-974f-4a9f-8e9c-e83001819942/iaudioclient-initialize-failure

		// initialization failed, try to use stereo waveformat
		m_wfx->nChannels = 2;
		m_wfx->nBlockAlign = (2 * m_wfx->wBitsPerSample) / 8;
		m_wfx->nAvgBytesPerSec = m_wfx->nSamplesPerSec * m_wfx->nBlockAlign;

		hr = m_clAudio->Initialize(AUDCLNT_SHAREMODE_SHARED, m_port == PORT_OUTPUT ? AUDCLNT_STREAMFLAGS_LOOPBACK : 0, REFTIMES_PER_SEC/*hnsRequestedDuration*/, 0, m_wfx, NULL);
		if (hr != S_OK)	{
			// stereo waveformat didnt work either, throw an error
			std::cerr << "WARNING: Failed to initialize audio client.\n";
		}
	}
	EXIT_ON_ERROR(hr);

	// initialize the audio capture client
	hr = m_clAudio->GetService(IID_IAudioCaptureClient, (void**)&m_clCapture);
	if (hr != S_OK)	{
		std::cerr << "WARNING: Failed to create audio capture client.\n";
	}
	EXIT_ON_ERROR(hr);

	// start the stream
	hr = m_clAudio->Start();
	if (hr != S_OK)	{
		std::cerr << "WARNING: Failed to start the stream.\n";
	}
	EXIT_ON_ERROR(hr);

	// initialize the watchdog timer
	QueryPerformanceCounter(&m_pcFill);

	return S_OK;

Exit:
	DeviceRelease();
	return hr;
}


/**
 * Release handles to audio resources.  (except the enumerator)
 */
void Measure::DeviceRelease() {
#if (WINDOWS_BUG_WORKAROUND)
	if (m_clBugAudio) {
		std::cout << "DEBUG: Releasing dummy stream audio device.\n";
		m_clBugAudio->Stop();
	}
	SAFE_RELEASE(m_clBugRender);
	SAFE_RELEASE(m_clBugAudio);
#endif

	if (m_clAudio) {
		std::cout << "DEBUG: Releasing audio device.\n";
		m_clAudio->Stop();
	}

	SAFE_RELEASE(m_clCapture);

	if (m_wfx) {
		CoTaskMemFree(m_wfx);
		m_wfx = NULL;
	}

	SAFE_RELEASE(m_clAudio);
	SAFE_RELEASE(m_dev);

	//for (int iChan = 0; iChan < Measure::MAX_CHANNELS; ++iChan) {
		if (m_fftCfg/*[iChan]*/) kiss_fftr_free(m_fftCfg/*[iChan]*/);
		m_fftCfg/*[iChan]*/ = NULL;

		if (m_fftIn/*[iChan]*/) free(m_fftIn/*[iChan]*/);
		m_fftIn/*[iChan]*/ = NULL;

		if (m_fftOut/*[iChan]*/) free(m_fftOut/*[iChan]*/);
		m_fftOut/*[iChan]*/ = NULL;

		//if (m_bandOut/*[iChan]*/) free(m_bandOut/*[iChan]*/);
		//m_bandOut/*[iChan]*/ = NULL;
	//}

	/*if (m_bandFreq)	{
		free(m_bandFreq);
		m_bandFreq = NULL;
	}*/

	if (m_fftTmpOut) {
		free(m_fftTmpOut);
		free(m_fftTmpIn);
		free(m_fftKWdw);
		m_fftTmpOut = NULL;
		m_fftTmpIn = NULL;
		m_fftKWdw = NULL;
		kiss_fft_cleanup();
	}

	m_devName[0] = '\0';
	m_format = FMT_INVALID;
}