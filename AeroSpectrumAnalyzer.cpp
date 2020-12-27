// Run program: Ctrl + F5 or Debug > Start Without Debugging menu
// Debug program: F5 or Debug > Start Debugging menu

// Tips for Getting Started: 
//   1. Use the Solution Explorer window to add/manage files
//   2. Use the Team Explorer window to connect to source control
//   3. Use the Output window to see build output and other messages
//   4. Use the Error List window to view errors
//   5. Go to Project > Add New Item to create new code files, or Project > Add Existing Item to add existing code files to the project
//   6. In the future, to open this project again, go to File > Open > Project and select the .sln file

#include <iostream>
#include <Windows.h>
//#include <cstdio>
#include <AudioClient.h>
#include <AudioPolicy.h>
#include <MMDeviceApi.h>
#include <FunctionDiscoveryKeys_devpkey.h>

#include <cmath>
#include <cassert>
#include <vector>

#include "kiss_fft130/kiss_fftr.h"

//#include "options.h"


// Overview: Audio level measurement from the Window Core Audio API
// See: http://msdn.microsoft.com/en-us/library/windows/desktop/dd370800%28v=vs.85%29.aspx

// REFERENCE_TIME time units per second and per millisecond
#define WINDOWS_BUG_WORKAROUND	0
#define REFTIMES_PER_SEC		10000000
#define TWOPI					(2*3.14159265358979323846)
#define EXIT_ON_ERROR(hres)		if(FAILED(hres)) { goto Exit; }
#define SAFE_RELEASE(p)			if((p)!=NULL) { (p)->Release(); (p) = NULL; }
#define CLAMP01(x) max(0.0, min(1.0, (x)))

#define EMPTY_TIMEOUT			0.500
#define DEVICE_TIMEOUT			1.500
#define QUERY_TIMEOUT			(1.0/60)
#define BANDS 19

enum Port {
	PORT_OUTPUT,
	PORT_INPUT,
};
/*enum Channel {
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
	TYPE_FFT, <----- ALWAYS THIS
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
};*/
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

//1:Input, 2:Output
Port					m_port = PORT_OUTPUT;						// port specifier (parsed from options)
//0:Left, 1:Right, Surround Options - 2:Center, 3:Sub, 4:Back Left, 5:Back Right, 6:Sub Left, 7:Sub Right, Other - 8:Sum, 9:Average
//Channel					m_channel = LEFT;					// channel specifier (parsed from options)
//Type					m_type = TYPE_FFT;						// data type specifier (parsed from options)
Format					m_format = FMT_PCM_F32;					// format specifier (detected in init)
//int						m_envRMS[2];				// RMS attack/decay times in ms (parsed from options)
//int						m_envPeak[2];				// peak attack/decay times in ms (parsed from options)
int						m_envFFT[2] = { 300, 300 };				// FFT attack/decay times in ms (parsed from options), 15, 250
int						m_fftSize = 1024;					// size of FFT (parsed from options), powers of 2 work best
int						m_fftOverlap = 0;				// number of samples between FFT calculations, between 0 and fftSize
//int						m_nBands = BANDS;					// number of frequency bands (parsed from options)
//double					m_gainRMS;					// RMS gain (parsed from options)
//double					m_gainPeak;					// peak gain (parsed from options)
double					m_freqMin = 20;					// min freq for band measurement
double					m_freqMax = 20000;					// max freq for band measurement
double					m_sensitivity = 35.0;				// dB range for FFT/Band return values (parsed from options), min 1

IMMDeviceEnumerator* m_enum;						// audio endpoint enumerator
IMMDevice* m_dev;						// audio endpoint device
WAVEFORMATEX* m_wfx;						// audio format info
IAudioClient* m_clAudio;					// audio client instance
IAudioCaptureClient* m_clCapture;				// capture client instance

const CLSID CLSID_MMDeviceEnumerator = __uuidof(MMDeviceEnumerator);
const IID IID_IMMDeviceEnumerator = __uuidof(IMMDeviceEnumerator);
const IID IID_IAudioClient = __uuidof(IAudioClient);
const IID IID_IAudioCaptureClient = __uuidof(IAudioCaptureClient);
const IID IID_IAudioRenderClient = __uuidof(IAudioRenderClient);

#if (WINDOWS_BUG_WORKAROUND)
IAudioClient* m_clBugAudio;				// audio client for dummy silent channel
IAudioRenderClient* m_clBugRender;				// render client for dummy silent channel
#endif
//WCHAR					m_reqID[64];				// requested device ID (parsed from options)
WCHAR					m_devName[64];				// device friendly name (detected in init)
//float					m_kRMS[2];					// RMS attack/decay filter constants
//float					m_kPeak[2];					// peak attack/decay filter constants
float					m_kFFT[2];					// FFT attack/decay filter constants
double					m_pcMult;					// performance counter inv frequency
LARGE_INTEGER			m_pcFill;					// performance counter on last full buffer
LARGE_INTEGER			m_pcPoll;					// performance counter on last device poll
kiss_fftr_cfg			m_fftCfg;// [MAX_CHANNELS] ;		// FFT states for each channel
float* m_fftIn;// [MAX_CHANNELS] ;		// buffer for each channel's FFT input
float* m_fftOut;// [MAX_CHANNELS] ;		// buffer for each channel's FFT output
float* m_fftKWdw;					// window function coefficients
float* m_fftTmpIn;					// temp FFT processing buffer
kiss_fft_cpx* m_fftTmpOut;				// temp FFT processing buffer
int	m_fftBufW;					// write index for input ring buffers
int	m_fftBufP;					// decremental counter - process FFT at zero
float m_bandFreq[BANDS]; // buffer of band max frequencies
float m_bandOut[BANDS]; // buffer of band values

void DeviceRelease() {
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

	if (m_fftCfg) kiss_fftr_free(m_fftCfg);
	m_fftCfg = NULL;

	if (m_fftIn) free(m_fftIn);
	m_fftIn = NULL;

	if (m_fftOut) free(m_fftOut);
	m_fftOut = NULL;

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

HRESULT	DeviceInit() {
	HRESULT hr;

	// get the device handle
	assert(m_enum && !m_dev);

	// get the default ID
	hr = m_enum->GetDefaultAudioEndpoint(m_port == PORT_OUTPUT ? eRender : eCapture, eConsole, &m_dev);
	
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
	if (hr != S_OK) {
		std::cerr << "WARNING: Failed to create audio client for Windows bug workaround.\n";
	}
#endif

	// get the main audio client
	hr = m_dev->Activate(IID_IAudioClient, CLSCTX_ALL, NULL, (void**)&m_clAudio);
	if (hr != S_OK) {
		std::cerr << "WARNING: Failed to create audio client.\n";
	}

	EXIT_ON_ERROR(hr);

	// parse audio format - Note: not all formats are supported.
	hr = m_clAudio->GetMixFormat(&m_wfx);
	EXIT_ON_ERROR(hr);

	switch (m_wfx->wFormatTag) {
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
	m_fftCfg = kiss_fftr_alloc(m_fftSize, 0, NULL, NULL);
	m_fftIn = (float*)calloc(m_fftSize * sizeof(float), 1);
	m_fftOut = (float*)calloc(m_fftSize * sizeof(float), 1);
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
	if (hr != S_OK) {
		// Compatibility with the Nahimic audio driver
		// https://github.com/rainmeter/rainmeter/commit/0a3dfa35357270512ec4a3c722674b67bff541d6
		// https://social.msdn.microsoft.com/Forums/windowsdesktop/en-US/bd8cd9f2-974f-4a9f-8e9c-e83001819942/iaudioclient-initialize-failure

		// initialization failed, try to use stereo waveformat
		m_wfx->nChannels = 2;
		m_wfx->nBlockAlign = (2 * m_wfx->wBitsPerSample) / 8;
		m_wfx->nAvgBytesPerSec = m_wfx->nSamplesPerSec * m_wfx->nBlockAlign;

		hr = m_clAudio->Initialize(AUDCLNT_SHAREMODE_SHARED, m_port == PORT_OUTPUT ? AUDCLNT_STREAMFLAGS_LOOPBACK : 0, REFTIMES_PER_SEC/*hnsRequestedDuration*/, 0, m_wfx, NULL);
		if (hr != S_OK) {
			// stereo waveformat didnt work either, throw an error
			std::cerr << "WARNING: Failed to initialize audio client.\n";
		}
	}
	EXIT_ON_ERROR(hr);

	// initialize the audio capture client
	hr = m_clAudio->GetService(IID_IAudioCaptureClient, (void**)&m_clCapture);
	if (hr != S_OK) {
		std::cerr << "WARNING: Failed to create audio capture client.\n";
	}
	EXIT_ON_ERROR(hr);

	// start the stream
	hr = m_clAudio->Start();
	if (hr != S_OK) {
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


void Initialize() {
	//m_reqID[0] = '\0';
	m_devName[0] = '\0';
	for (int i = 0; i < BANDS; i++) {
		m_bandOut[i] = 0;
		m_bandFreq[i] = 0;
	}

	// create the enumerator
	if (CoCreateInstance(CLSID_MMDeviceEnumerator, NULL, CLSCTX_ALL, IID_IMMDeviceEnumerator, (void**)&m_enum) == S_OK) {
		// init the device (ok if it fails - it'll keep checking during Update)
		DeviceInit();
		return;
	}

	// regenerate filter constants
	if (m_wfx) {
		const double freq = m_wfx->nSamplesPerSec;
		//		m_kRMS[0] = (float)exp(log10(0.01) / (freq * (double)m_envRMS[0] * 0.001));
		//		m_kRMS[1] = (float)exp(log10(0.01) / (freq * (double)m_envRMS[1] * 0.001));
		//		m_kPeak[0] = (float)exp(log10(0.01) / (freq * (double)m_envPeak[0] * 0.001));
		//		m_kPeak[1] = (float)exp(log10(0.01) / (freq * (double)m_envPeak[1] * 0.001));
		m_kFFT[0] = (float)exp(log10(0.01) / (freq / (m_fftSize - m_fftOverlap) * (double)m_envFFT[0] * 0.001));
		m_kFFT[1] = (float)exp(log10(0.01) / (freq / (m_fftSize - m_fftOverlap) * (double)m_envFFT[1] * 0.001));
	}

	SAFE_RELEASE(m_enum)
}

void Update() {
	LARGE_INTEGER	pcCur;
	QueryPerformanceCounter(&pcCur);

	// query the buffer
	if (m_clCapture && ((pcCur.QuadPart - m_pcPoll.QuadPart) * m_pcMult) >= QUERY_TIMEOUT) {
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
				(m_fftIn)[m_fftBufW] = (m_format == FMT_PCM_F32) ? *sF32++ : ((float)*sI16++ * 1.0f / 0x7fff);
				m_fftBufW = (m_fftBufW + 1) % m_fftSize;
				// if overlap limit reached, process FFTs for each channel
				if (!--m_fftBufP) {
					if (!(flags & AUDCLNT_BUFFERFLAGS_SILENT)) {
						// copy from the ring buffer to temp space
						memcpy(&m_fftTmpIn[0], &(m_fftIn)[m_fftBufW], (m_fftSize - m_fftBufW) * sizeof(float));
						memcpy(&m_fftTmpIn[m_fftSize - m_fftBufW], &m_fftIn[0], m_fftBufW * sizeof(float));
						// apply the windowing function
						for (int iBin = 0; iBin < m_fftSize; ++iBin) {
							m_fftTmpIn[iBin] *= m_fftKWdw[iBin];
						}
						kiss_fftr(m_fftCfg, m_fftTmpIn, m_fftTmpOut);
					}
					else {
						memset(m_fftTmpOut, 0, m_fftSize * sizeof(kiss_fft_cpx));
					}
					// filter the bin levels as with peak measurements
					for (int iBin = 0; iBin < m_fftSize; ++iBin) {
						float x0 = (m_fftOut)[iBin];
						float x1 = (m_fftTmpOut[iBin].r * m_fftTmpOut[iBin].r + m_fftTmpOut[iBin].i * m_fftTmpOut[iBin].i) * scalar;
						x0 = x1 + m_kFFT[(x1 < x0)] * (x0 - x1);
						(m_fftOut)[iBin] = x0;
					}
					m_fftBufP = m_fftSize - m_fftOverlap;
				}

				// integrate FFT results into log-scale frequency bands
				const float		df = (float)m_wfx->nSamplesPerSec / m_fftSize;
				const float		scalar = 2.0f / (float)m_wfx->nSamplesPerSec;
				memset(m_bandOut, 0, BANDS * sizeof(float));
				int			iBin = 0;
				int			iBand = 0;
				float		f0 = 0.0f;
				while (iBin <= (m_fftSize / 2) && iBand < BANDS) {
					float	fLin1 = ((float)iBin + 0.5f) * df;
					float	fLog1 = m_bandFreq[iBand];
					float	x = (m_fftOut)[iBin];
					float& y = (m_bandOut)[iBand];
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
	}
	else if (!m_clCapture && (((pcCur.QuadPart - m_pcPoll.QuadPart) * m_pcMult) >= DEVICE_TIMEOUT)) {
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
	}*/
}


int main() {
	CoInitialize(nullptr); // NULL if using older VC++
	Initialize();
	Update();

	std::cout << m_devName;

	/*while (1) {
		Update();
		for (int i = 0; i < BANDS; i++) {
			float& y = m_bandOut[i];
			std::cout << y << ",";
		}
		std::cout << std::endl;
	}*/

	return 0;
}









/*#include "measure.h"

int main() {
	Options* o = new Options("Options.txt");
	Measure* m = new Measure(o);
	
	while (1) {
		m->Update();
		for (int i = 0; i < BANDS; i++) {
			float& y = m->m_bandOut[i];
			std::cout << y << ",";
		}
		std::cout << std::endl;
	}

    return 0;
}*/