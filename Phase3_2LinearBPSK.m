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
encNumberOfBits = numberOfBits * 7 / 4;

randomInputBits = randi([0, 1], [1, numberOfBits]);
%%encRandomInputBits = randi([0, 1], [1, encNumberOfBits]);

pol = cyclpoly(7,4);
parmat = cyclgen(7,pol);
genmat = gen2par(parmat);
encRandomInputBits = encode(randomInputBits,7,4,'linear/binary',genmat);

% Get all the timepoints we used to sample the data, it starts at half of the
% sampling intervel (mid point of somallest definition)
timePoint = 1/(2 * sampleFreq):  1/sampleFreq : numberOfBits/dataRate;
encTimePoint = 1/(2 * sampleFreq):  1/sampleFreq : encNumberOfBits/dataRate;


% Carrier singal with marked timepoints
carrierSingal = cos(2 * pi * carrierFreq * timePoint);
encCarrierSingal =  cos(2 * pi * carrierFreq * encTimePoint);
% obtained the sampled signal by extending the orignal signal
sampledSingal = kron(randomInputBits, ones(1, sampleFreq/dataRate));
encSampledSingal = kron(encRandomInputBits, ones(1, sampleFreq/dataRate));


% Sampled OOK signal, 1 as 1 and 0 remains as 0
sampledOOK = sampledSingal .* carrierSingal;
encSampledOOK = encSampledSingal .* encCarrierSingal;
% Sampled BPSK signal needs to change 0 to -1 as phase changes 180 degress
% while 1 remains as 1
sampledBPSK = (2 * sampledSingal - 1) .* carrierSingal;
encSampledBPSK = (2 * encSampledSingal - 1) .* encCarrierSingal;



%---------------------------------------------------------------------
%Testing Code
%---------------------------------------------------------------------

% %Create Output Results Array
% bitErrorRateOutput = zeros(1,11); %Y axis
% encBitErrorRateOutput = zeros(1,11); %Y axis
% SNRAxis = zeros(1,11); %X axis
% %Carry out test considering different SNR values from 0 dB to 50 dB (in multiples of 5 dB)
% for i = 0:10
    % SNR = 5 * i;
    % numberOfTrails = 50; % For every SNR we transmit 100 times and take average to get errorrate
    % bitErrorRateTrails = zeros(numberOfTrails, 1);
	% encBitErrorRateTrails = zeros(numberOfTrails, 1);
    % SNRAxis(i+1) = SNR;
    % for j = 1:numberOfTrails
        
        % decodedSignal = tansmitAndDemodulateFunction(sampledOOK, signalVariance, SNR, 0.5, numberOfBits, sampleFreq, dataRate, carrierSingal,b,a);
        % bitErrorRateTrails(j) = calculateBitError(randomInputBits, decodedSignal);
		
		% encTransmittedSignal = tansmitAndDemodulateFunction(encSampledOOK, signalVariance, SNR, 0.5, encNumberOfBits, sampleFreq, dataRate, encCarrierSingal,b,a);
        % encDecodedSignal = decode(encTransmittedSignal,7,4,'linear/binary',genmat);
		% encBitErrorRateTrails(j) = calculateBitError(randomInputBits, encDecodedSignal);
        
    % end
    % %Calculate the mean errorrate as the final result of certain SNR
    % bitErrorRateOutput(i+1) = mean(bitErrorRateTrails);
	% encBitErrorRateOutput(i+1) = mean(encBitErrorRateTrails);
% end

% semilogy(SNRAxis, bitErrorRateOutput);
% hold on
% semilogy(SNRAxis, encBitErrorRateOutput);
% ylim([10^(-5) 10^1]);
% xlim([0 50]);
% hold on


%Create Output Results Array
bitErrorRateOutputBPSK = zeros(1,11); %Y axis
encBitErrorRateOutputBPSK = zeros(1,11); %Y axis
SNRAxis = zeros(1,11); %X axis
%Carry out test considering different SNR values from 0 dB to 50 dB (in multiples of 5 dB)
for i = 0:10
    SNR = 5 * i;
    numberOfTrails = 50; % For every SNR we transmit 50  times and take average to get errorrate
    bitErrorRateTrails = zeros(numberOfTrails, 1);
	encBitErrorRateTrails = zeros(numberOfTrails, 1);
    SNRAxis(i+1) = SNR;
    for j = 1:numberOfTrails
        
        decodedSignal = tansmitAndDemodulateFunction(sampledBPSK, signalVariance, SNR, 0, numberOfBits, sampleFreq, dataRate, carrierSingal,b,a);
        bitErrorRateTrails(j) = calculateBitError(randomInputBits, decodedSignal);
        
		encTransmittedSignal = tansmitAndDemodulateFunction(encSampledBPSK, signalVariance, SNR, 0, encNumberOfBits, sampleFreq, dataRate, encCarrierSingal,b,a);
        encDecodedSignal = decode(encTransmittedSignal,7,4,'linear/binary',genmat);
		encBitErrorRateTrails(j) = calculateBitError(randomInputBits, encDecodedSignal);
        
    end
    %Calculate the mean errorrate as the final result of certain SNR
    bitErrorRateOutputBPSK(i+1) = mean(bitErrorRateTrails); 
	encBitErrorRateOutputBPSK(i+1) = mean(encBitErrorRateTrails);
end


semilogy(SNRAxis, bitErrorRateOutputBPSK);
hold on
semilogy(SNRAxis, encBitErrorRateOutputBPSK);
hold on
ylim([0 1]);
xlim([0 50]);
hold on
ylabel('Log 10 Bit Error Rate') ;
hold on
title("Bit Error Rate vs SNR for BPSK with and without Linear coding");
legend({'y = AverageBPSK','y= AverageBPSK'},'Location','southeast')
xlabel('E_{b}/N_{0}') ;
ylabel('P_{e}') ;





%---------------------------------------------------------------------
%Function Code
%---------------------------------------------------------------------

function decodedSignal = tansmitAndDemodulateFunction(input, inputVariance, SNR, threshHold, numberOfBits, sampleFreq, dataRate, carrierSingal, filterB, filterA)
    % Generate equal number of noise samples
    % The generated noise should have normal distribution with zero mean and unit variance (use the function randn in MATLAB).
    mean = 0;
    
    
    % Change the noise variance with respect to SNR (signal to noise ratio) value
    % Use SNR value to generate noise variance
    % SNR (in dB) = 10log10 (S/N) where S is the Signal power (or variance) and N is the Noise power (or variance)
    noiseVariance = inputVariance / 10^(SNR/10);
    noiseSTD = sqrt(noiseVariance);
    noiseBits = noiseSTD.*randn( 1, numberOfBits * sampleFreq/dataRate) + mean;
    
    %Find output with noise data
    outputWithNoise = input + noiseBits;
    
    %DemoduleteSignal and pass it will low-pass-filter
    demodulatedSignal = outputWithNoise.*(2 * carrierSingal);
    filteredSignal = filtfilt(filterB, filterA, demodulatedSignal);
    decodedSignal = zeros(1,numberOfBits);
    
    for i = 1:1:numberOfBits
        decodedSignalSample = filteredSignal(1 / 2 * sampleFreq/dataRate + (i - 1) * sampleFreq/dataRate);
        if decodedSignalSample > threshHold
            decodedSignal(i) = 1;
        else
            decodedSignal(i) = 0;        
        end
    end
    
    
    
    
    
end


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