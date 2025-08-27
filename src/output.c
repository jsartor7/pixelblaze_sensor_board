
#include "main.h"

//an fft of 512 gives us 256 buckets of frequency info. at 20khz, each has ~38Hz
//we don't really want all 256 buckets of frequency info
//if we compressed this linearly, we'd lose a lot of low/mid tone info
//also, in order to capture low frequency stuff below 38Hz would require a much larger fft
//so we can combine a downsampled 400hz using a smaller fft for low frequency stuff with the 20khz stuff

//the first 6 buckets are dedicated to low frequency audio from 12.5-162.5 Hz




// 16.16 fixed-point arithmetic macros
#define FIX16_ONE (1 << 16)
#define FIX16_FROM_INT(x) ((x) << 16)
#define FIX16_TO_INT(x) ((x) >> 16)
#define FIX16_MUL(x, y) (((int64_t)(x) * (y)) >> 16)
#define FIX16_DIV(x, y) (((int64_t)(x) << 16) / (y))
#define FIX16_FROM_FLOAT(x) ((int32_t)((x) * 65536.0f))

const int32_t numFreqs = 88;

//these are all of the notes on a piano, the first in each line (except for the first) is an C
const int32_t noteFrequencies[88] = {
	FIX16_FROM_FLOAT(27.50), FIX16_FROM_FLOAT(29.14), FIX16_FROM_FLOAT(30.87),
	FIX16_FROM_FLOAT(32.70), FIX16_FROM_FLOAT(34.65), FIX16_FROM_FLOAT(36.71), FIX16_FROM_FLOAT(38.89), FIX16_FROM_FLOAT(41.20), FIX16_FROM_FLOAT(43.65), FIX16_FROM_FLOAT(46.25), FIX16_FROM_FLOAT(49), FIX16_FROM_FLOAT(51.91), FIX16_FROM_FLOAT(55), FIX16_FROM_FLOAT(58.27), FIX16_FROM_FLOAT(61.74),
	FIX16_FROM_FLOAT(65.41), FIX16_FROM_FLOAT(69.30), FIX16_FROM_FLOAT(73.42), FIX16_FROM_FLOAT(77.78), FIX16_FROM_FLOAT(82.41), FIX16_FROM_FLOAT(87.31), FIX16_FROM_FLOAT(92.50), FIX16_FROM_FLOAT(98), FIX16_FROM_FLOAT(103.83), FIX16_FROM_FLOAT(110), FIX16_FROM_FLOAT(116.54), FIX16_FROM_FLOAT(123.47),
	FIX16_FROM_FLOAT(130.81), FIX16_FROM_FLOAT(138.59), FIX16_FROM_FLOAT(146.83), FIX16_FROM_FLOAT(155.56), FIX16_FROM_FLOAT(164.81), FIX16_FROM_FLOAT(174.61), FIX16_FROM_FLOAT(185), FIX16_FROM_FLOAT(196), FIX16_FROM_FLOAT(207.65), FIX16_FROM_FLOAT(220), FIX16_FROM_FLOAT(233.08), FIX16_FROM_FLOAT(246.94),
	FIX16_FROM_FLOAT(261.63), FIX16_FROM_FLOAT(277.18), FIX16_FROM_FLOAT(293.66), FIX16_FROM_FLOAT(311.13), FIX16_FROM_FLOAT(329.63), FIX16_FROM_FLOAT(349.23), FIX16_FROM_FLOAT(369.99), FIX16_FROM_FLOAT(392), FIX16_FROM_FLOAT(415.30), FIX16_FROM_FLOAT(440), FIX16_FROM_FLOAT(466.16), FIX16_FROM_FLOAT(493.88),
	FIX16_FROM_FLOAT(523.25), FIX16_FROM_FLOAT(554.37), FIX16_FROM_FLOAT(587.33), FIX16_FROM_FLOAT(622.25), FIX16_FROM_FLOAT(659.25), FIX16_FROM_FLOAT(698.46), FIX16_FROM_FLOAT(739.99), FIX16_FROM_FLOAT(783.99), FIX16_FROM_FLOAT(830.61), FIX16_FROM_FLOAT(880), FIX16_FROM_FLOAT(932.33), FIX16_FROM_FLOAT(987.77),
	FIX16_FROM_FLOAT(1046.50), FIX16_FROM_FLOAT(1108.73), FIX16_FROM_FLOAT(1174.66), FIX16_FROM_FLOAT(1244.51), FIX16_FROM_FLOAT(1318.51), FIX16_FROM_FLOAT(1396.91), FIX16_FROM_FLOAT(1479.98), FIX16_FROM_FLOAT(1567.98), FIX16_FROM_FLOAT(1661.22), FIX16_FROM_FLOAT(1760), FIX16_FROM_FLOAT(1864.66), FIX16_FROM_FLOAT(1975.53),
	FIX16_FROM_FLOAT(2093.00), FIX16_FROM_FLOAT(2217.46), FIX16_FROM_FLOAT(2349.32), FIX16_FROM_FLOAT(2489.02), FIX16_FROM_FLOAT(2637.02), FIX16_FROM_FLOAT(2793.83), FIX16_FROM_FLOAT(2959.96), FIX16_FROM_FLOAT(3135.96), FIX16_FROM_FLOAT(3322.44), FIX16_FROM_FLOAT(3520), FIX16_FROM_FLOAT(3729.31), FIX16_FROM_FLOAT(3951.07),
	FIX16_FROM_FLOAT(4186.01)
};

//this function takes a value at a frequency and converts it to a note and then adds it to our note values
void updateNoteMags(uint16_t * magnitude, int downsampleFreq, uint16_t numSteps, uint16_t * noteMags) {

	for (int i = 1; i < numSteps/2; i++)
	{
		int32_t stepSize = FIX16_DIV(FIX16_FROM_INT(downsampleFreq), FIX16_FROM_INT(numSteps));
		int32_t currFreq = FIX16_MUL(stepSize, FIX16_FROM_INT(i));

		//find the index j where our fourier frequency is greater than a note frequency
		int j = 0;
		while(j < (int) numFreqs && noteFrequencies[j] < currFreq)
		{
			j++;
			//noteMags[j%12] += j*1000;
		}
		//after this loop, j is the first freq index that is greater than currFreq

		uint8_t noteNum = 12;
		int32_t distance = FIX16_FROM_INT(100);
		int32_t reject_distance = FIX16_FROM_INT(100);
		//ensure we're not just above or below the whole range
		if(j > 0 && j < (numFreqs-2) && currFreq > (noteFrequencies[0] - stepSize))
		{
			int32_t distanceUp = (noteFrequencies[j] - currFreq);
			int32_t distanceDown = (currFreq - noteFrequencies[j-1]);

			//decide which one is closer
			if(distanceUp < 0 || distanceDown < 0)
			{
				noteNum = 12;
			}
			else if(distanceUp < distanceDown)
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
		if( distance < stepSize/2 && distance < reject_distance/2){

			if (currFreq > FIX16_FROM_INT(1600))
			{
				noteMags[noteNum] += magnitude[i] << 3;
			}
			else if (currFreq > FIX16_FROM_INT(800))
			{
				noteMags[noteNum] += magnitude[i] << 2;
			}
			else if (currFreq > FIX16_FROM_INT(400))
			{
				noteMags[noteNum] += magnitude[i] << 1;
			}
			else
			{
				noteMags[noteNum] += magnitude[i];
			}
		}
	}

}




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


//void processSensorData(int16_t * audioBuffer, int16_t * audioLowHzBuffer, int16_t * audioMidHzBuffer, volatile uint16_t adcBuffer[7], volatile int16_t accelerometer[3]) {
void processSensorData(int16_t * audioBuffer, int16_t * audioMidHzBuffer, volatile uint16_t adcBuffer[7], volatile int16_t accelerometer[3]) {
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
	//fftRealWindowedMagnitude(audioLowHzBuffer, &magnitude[0], LOW_NLOG2, &lowEnergy);

	//updateNoteMags(&magnitude[0], 40000 / LOW_N_DOWNSAMPLE, LOW_N, &noteMags[0]);


	//if()
	fftRealWindowedMagnitude(audioMidHzBuffer, &magnitude[0], MID_NLOG2, &lowEnergy);

	updateNoteMags(&magnitude[0], 40000 / MID_N_DOWNSAMPLE, MID_N, &noteMags[0]);


	for (int i = 0; i < 32; i++)
	{
		if(i < 12)
		{
		WRITEOUT(noteMags[i%12]);
		}
		else
		{
			WRITEOUT(noteMags[0]);
			//WRITEOUT(magnitude[(i-12)*MID_N])
		}
		//WRITEOUT(magnitude[i]);
	}

/*
	for (int i = 0, k = lowFrequencyMap[0]; i < 6; i++) {
		int top = lowFrequencyMap[i] + 1;
		uint16_t max = 0;
		for (; k < top; k++) {
			max = magnitude[k] > max ? magnitude[k] : max;
		}
		WRITEOUT(max);
	}
*/
	//do high frequency stuff
	//fftRealWindowedMagnitude(audioBuffer, &magnitude[0], HIGH_NLOG2, &energyAverage);

	//run through and get maxFrequency info
	for (int i = 1; i < HIGH_N/2; i++) {
		if (magnitude[i] > maxFrequencyMagnitude) {
			maxFrequencyMagnitude = magnitude[i];
			maxFrequencyIndex = i;
		}
	}
/*
	//write out high frequency stuff
	for (int i = 0, k = highFrequencyMap[0]; i < 26; i++) {
		int top = highFrequencyMap[i] + 1;
		uint16_t max = 0;
		for (; k < top; k++) {
			max = magnitude[k] > max ? magnitude[k] : max;
		}
		WRITEOUT(max);
	}
	*/

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
