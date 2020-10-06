%---------------------------------------------------------------------
%Preparation Code
%---------------------------------------------------------------------

% Define the carrier frequency as 10 KHz
carrierFreq = 100000;
% Define the carrier signal as 16 times oversampled 
sampleFreq = carrierFreq * 16;
% Define the dataRate of original signal 1kbps
dataRate = 1000;
% Define low pass butter filter, 6th order with cut of freq of 0.2
[b, a] = butter(6, 0.2);
% Define amp
amplitude = 1;
% Define signal (the input data) has unit power. That is, S=1
signalVariance = 1;


% Generate random binary digits (0 or 1)for 1024 bits
numberOfBits = 1024;
randomInputBits = randi([0, 1], [1, numberOfBits]);


% Get all the timepoints we used to sample the data, it starts at half of the
% sampling intervel (mid point of somallest definition)
timePoint = 1/(2 * sampleFreq):  1/sampleFreq : numberOfBits/dataRate;


% Carrier singal with marked timepoints
carrierSingal = cos(2 * pi * carrierFreq * timePoint);
% obtained the sampled signal by extending the orignal signal
sampledSingal = kron(randomInputBits, ones(1, sampleFreq/dataRate));

% Sampled BPSK signal needs to change 0 to -1 as phase changes 180 degress
% while 1 remains as 1
sampledBPSK = (2 * sampledSingal - 1) .* carrierSingal;

%---------------------------------------------------------------------
%Testing Code
%---------------------------------------------------------------------

%Create Output Results Array
bitErrorRateOutput = zeros(1,11); %Y axis
SNRAxis = zeros(1,11); %X axis
%Carry out test considering different SNR values from 0 dB to 50 dB (in multiples of 5 dB)
SNR = 8;
    
% Generate equal number of noise samples
% The generated noise should have normal distribution with zero mean and unit variance (use the function randn in MATLAB).
mean = 0;
    
% Change the noise variance with respect to SNR (signal to noise ratio) value
% Use SNR value to generate noise variance
% SNR (in dB) = 10log10 (S/N) where S is the Signal power (or variance) and N is the Noise power (or variance)
noiseVariance = signalVariance / 10^(SNR/10);
noiseSTD = sqrt(noiseVariance);
noiseBits = noiseSTD.*randn( 1, numberOfBits * sampleFreq/dataRate) + mean;

%Find output with noise data
outputWithNoise = sampledOOK + noiseBits;

%DemoduleteSignal and pass it will low-pass-filter
demodulatedSignal = outputWithNoise.*(2 * carrierSingal);
filteredSignal = filtfilt(b, a, demodulatedSignal);
decodedSignal = zeros(1,numberOfBits);


for i = 1:1:numberOfBits
    decodedSignalSample = filteredSignal(1 / 2 * sampleFreq/dataRate + (i - 1) * sampleFreq/dataRate);
    if decodedSignalSample > 0
        decodedSignal(i) = 1;
    else
        decodedSignal(i) = 0;        
    end
end




bitErrorRate = calculateBitError(randomInputBits, decodedSignal);
decodedSignal = kron(decodedSignal, ones(1, sampleFreq/dataRate));


%---------------------------------------------------------------------
%Plotting Code
%---------------------------------------------------------------------

ts1 = timeseries(sampledSingal,timePoint);
ts1.Name = 'Sampled Signal';
subplot(5, 1, 1);
plot(ts1);
xlim([0 0.2]);
ylim([-2 2]);


ts2 = timeseries(sampledBPSK,timePoint);
ts2.Name = 'Modulated Signal';
subplot(5, 1, 2);
plot(ts2);
xlim([0 0.2]);
ylim([-2 2]);

ts3 = timeseries(outputWithNoise,timePoint);
ts3.Name = 'Received Signal';
subplot(5, 1, 3);
plot(ts3);
xlim([0 0.2]);
ylim([-4 4]);

ts4 = timeseries(demodulatedSignal,timePoint);
ts4.Name = 'Demodulated signal';
subplot(5, 1, 4);
plot(ts4);
xlim([0 0.2]);
ylim([-4 4]);

ts5 = timeseries(decodedSignal,timePoint);
ts5.Name = 'Decoded signal';
subplot(5, 1, 5);
plot(ts5);
xlim([0 0.2]);
ylim([-2 2]);

function bitErrorRate = calculateBitError(input, output)

    bitError = 0;
    bitErrorRate = 0;
    %verify equality of the length of input and output
    if length(input) ~= length(output)
        return
    end
        
    %Compute p(e)   
    for i = 1 : length(input)
        if (input(i) ~= output(i))
            bitError = bitError + 1;
        end
        
    bitErrorRate= bitError/length(input);
    end
    
    
end