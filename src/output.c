
#include "main.h"

//an fft of 512 gives us 256 buckets of frequency info. at 20khz, each has ~38Hz
//we don't really want all 256 buckets of frequency info
//if we compressed this linearly, we'd lose a lot of low/mid tone info
//also, in order to capture low frequency stuff below 38Hz would require a much larger fft
//so we can combine a downsampled 400hz using a smaller fft for low frequency stuff with the 20khz stuff

//the first 6 buckets are dedicated to low frequency audio from 12.5-162.5 Hz

const float numFreqs = 88;

//first row is C, etc.
//27.50 is lowest A on a piano, 4186 is highest C
const float noteFrequencies[88] = {
		//16.35,	17.32,	18.35,	19.45,	20.60,	21.83,	23.12,	24.50,	25.96,
		27.50,	29.14,	30.87,
		32.70,	34.65,	36.71,	38.89,	41.20,	43.65,	46.25,	49,	51.91,	55,	58.27,	61.74,
		65.41,	69.30,	73.42,	77.78,	82.41,	87.31,	92.50,	98,	103.83,	110,	116.54,	123.47,
		130.81,	138.59,	146.83,	155.56,	164.81,	174.61,	185,	196,	207.65,	220,	233.08,	246.94,
		261.63,	277.18,	293.66,	311.13,	329.63,	349.23,	369.99,	392,	415.30,	440,	466.16,	493.88,
		523.25,	554.37,	587.33,	622.25,	659.25,	698.46,	739.99,	783.99,	830.61,	880,	932.33,	987.77,
		1046.50, 1108.73, 1174.66, 1244.51, 1318.51, 1396.91, 1479.98, 1567.98, 1661.22, 1760, 1864.66, 1975.53,
		2093.00, 2217.46, 2349.32, 2489.02, 2637.02, 2793.83, 2959.96, 3135.96, 3322.44, 3520, 3729.31, 3951.07,
		4186.01
};


//this is a list of indices which refer to places within some other list somewhere else.
//I think we should probably make the low frequency one bigger and the other one smaller.

const uint8_t lowFrequencyMap[6] = {
		3, 4, 6, 8, 10, 13
};
//the next 26 buckets are dedicated to higher frequency audio from 195Hz to just under the nyquist limit of 10khz
const uint8_t highFrequencyMap[26] = {
		5, 6, 8, 10, 12, 15, 18, 22, 25, 30, 35, 40, 46, 53, 61, 70, 80, 92, 105, 119, 136, 154, 175, 199, 225, 255,
};

char outBuffer[100];
int outBufferLen;

#define WRITEOUT(v) {memcpy(out, &v, sizeof(v)); out+= sizeof(v);}

/*
 * Takes a real input, applies Hann window, calculates energyAverage
 * out must be the the same size as in as it is used for the imaginary part
 * after returning, the first half is filled with bucket magnitudes
 * magnitude is multiplied by 16 and saturates at 16 bits
 * m = log2(n)
 */
void fftRealWindowedMagnitude(int16_t * in, uint16_t * out, int m, uint16_t * energyAverage) {
	int16_t * imag = (int16_t *) out; //borrow out for imaginary part
	int n = 1 << m;

	//init with zeros
	for (int i = 0; i < n; i++) {
		imag[i] = 0;
	}
	uint32_t energyTotal = 0;
	for (int i = 0; i < n; i++) {
		energyTotal += abs(in[i]);

		//apply the hann windowing function, borrowing Sinewave LUT from fix_fft
		//the positive portion of Sinewave ranges from index 0-512
		// (i * 512) / n  == (i * 512) >> m
		int si = (i * 512) >> m ;
		in[i] = (Sinewave[si] * in[i]) >> 16;
	}
	*energyAverage = energyTotal >> m;

	//run the FFT (runs in place, overwriting both in and imag)
	fix_fft(in, imag, m, 0);

	//calculate the magnitude and store in out for only the first half
	int halfN = n>>1;
	for (int i = 0; i < halfN; i++) {
		//using the fix16_sqrt gives us a bit more resolution as we get
		//8 bits more using this over an integer sqrt
		int32_t t = fix16_sqrt(imag[i] * imag[i] + in[i] * in[i]);

		//we can't keep all those extra bits, but 4 of 8 seems like a good value
		//as this only overloads a little and only for REALLY LOUD inputs
		t >>= 4;
		if (t > 0xffff)
			t = 0xffff;
		out[i] = t;
	}
}

void updateNoteMags(uint16_t * magnitude, float downsampleFreq, uint16_t numSteps, uint16_t * noteMags) {

	for (int i = 0; i < numSteps/2; i++)
	{
		float stepSize = (float) downsampleFreq / (float) numSteps;
		float currFreq = stepSize * (float) i;
		//find the index j where our fourier frequency is greater than a note frequency
		int j = 0;
		while(j < (int) numFreqs && noteFrequencies[j] < currFreq)
		{
			j++;
			//noteMags[j%12] += j*1000;
		}

		uint8_t noteNum = 12;
		float distance = 100;
		float reject_distance = 100;
		//ensure we're not just above or below the whole range
		if(j < (numFreqs-2) && currFreq > (noteFrequencies[0]-stepSize))
		{
			float distanceUp = (noteFrequencies[j]-currFreq);
			float distanceDown = j > 0 ? (currFreq-noteFrequencies[j-1]) : 100;

			//decide which one is closer
			if(j == 0 || distanceUp < distanceDown)
			{
				distance = distanceUp;
				reject_distance = distanceDown;
				noteNum = j%12;
			}
			else
			{
				distance = distanceDown;
				reject_distance = distanceUp;
				noteNum = (j-1)%12;
			}
		}

		//a note was assigned, and it is reasonably close to the frequency in question
		if( distance < stepSize && distance < (reject_distance/2) && noteNum < 12)
		{
			noteMags[noteNum] += magnitude[i] / 4;
		}
	}

}

void processSensorData(int16_t * audioBuffer, int16_t * audio400HzBuffer, int16_t * audioMidHzBuffer, volatile uint16_t adcBuffer[7], volatile int16_t accelerometer[3]) {
	uint16_t magnitude[HIGH_N]; //temp and output from the fft
	uint16_t noteMags[12]; // these will be our main outputs. initializes to 0
	uint16_t lowEnergy;
	uint16_t energyAverage;
	int maxFrequencyIndex = 0;
	uint16_t maxFrequencyMagnitude = 0;
	uint16_t maxFrequencyHz;
	char * out = outBuffer;

	//start making output buffer
	WRITEOUT("SB1.0");

	for (int i = 0; i < 12; i++)
	{
		noteMags[i] = 0;
	}

	//do the low frequency stuff
	fftRealWindowedMagnitude(audio400HzBuffer, &magnitude[0], LOW_NLOG2, &lowEnergy);

	//updateNoteMags(&magnitude[0], 400, LOW_N, &noteMags[0]);

	//fftRealWindowedMagnitude(audioMidHzBuffer, &magnitude[0], MID_NLOG2, &lowEnergy);

	//updateNoteMags(&magnitude[0], 4000, MID_N, &noteMags[0]);

	//fftRealWindowedMagnitude(audio400HzBuffer, &magnitude[0], LOW_NLOG2, &lowEnergy);


	for (int i = 0; i < 32; i++)
	{
		//WRITEOUT(noteMags[i%12]);
		//WRITEOUT(magnitude[i]);
	}


	for (int i = 0, k = lowFrequencyMap[0]; i < 6; i++) {
		int top = lowFrequencyMap[i] + 1;
		uint16_t max = 0;
		for (; k < top; k++) {
			max = magnitude[k] > max ? magnitude[k] : max;
		}
		WRITEOUT(max);
	}

	//do high frequency stuff
	fftRealWindowedMagnitude(audioBuffer, &magnitude[0], HIGH_NLOG2, &energyAverage);

	//run through and get maxFrequency info
	for (int i = 1; i < HIGH_N/2; i++) {
		if (magnitude[i] > maxFrequencyMagnitude) {
			maxFrequencyMagnitude = magnitude[i];
			maxFrequencyIndex = i;
		}
	}

	//write out high frequency stuff
	for (int i = 0, k = highFrequencyMap[0]; i < 26; i++) {
		int top = highFrequencyMap[i] + 1;
		uint16_t max = 0;
		for (; k < top; k++) {
			max = magnitude[k] > max ? magnitude[k] : max;
		}
		WRITEOUT(max);
	}

	WRITEOUT(energyAverage);
	WRITEOUT(maxFrequencyMagnitude);
	maxFrequencyHz = (20000 * (int32_t)maxFrequencyIndex) / HIGH_N; //or 39.0625 per bin
	WRITEOUT(maxFrequencyHz);

	for (int i = 0; i < 3; i++) {
		int16_t v = accelerometer[i];
		WRITEOUT(v);
	}

	//light
	uint16_t v = adcBuffer[1]<<4;
	WRITEOUT(v);

	//a0
	v = adcBuffer[2]<<4;
	WRITEOUT(v);
	//a1
	v = adcBuffer[6]<<4;
	WRITEOUT(v);
	//a2
	v = adcBuffer[5]<<4;
	WRITEOUT(v);
	//a3
	v = adcBuffer[4]<<4;
	WRITEOUT(v);
	//a4
	v = adcBuffer[3]<<4;
	WRITEOUT(v);


	WRITEOUT("END");

	outBufferLen = out - outBuffer;
	writeToUsart((uint8_t *) outBuffer, outBufferLen);
}
