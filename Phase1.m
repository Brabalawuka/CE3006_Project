%---------------------------------------------------------------------
%Preparation Code
%---------------------------------------------------------------------

% Generate random binary digits (0 or 1)for 1024 bits
numberOfBits = 1024;
randomInputBits = randi([0, 1], [numberOfBits, 1]);

% Convert binary digits to ±1 (means 1 to +1 and 0 to -1).
randomInputBits(randomInputBits == 0) = -1;
% Assume signal (the input data) has unit power. That is, S=1
signalVariance = 1;

%---------------------------------------------------------------------
%Testing Code
%---------------------------------------------------------------------

%Create Output Results Array
bitErrorRateOutput = zeros(1,11); %Y axis
SNRAxis = zeros(1,11); %X axis

%Carry out test considering different SNR values from 0 dB to 50 dB (in multiples of 5 dB)
for i = 0:10
    SNR = 5 * i;
    numberOfTrails = 100; % For every SNR we transmit 100 times and take average to get errorrate
    bitErrorRateTrails = zeros(numberOfTrails, 1);
    SNRAxis(i+1) = SNR;
    for j = 1:numberOfTrails
        
        outputWithNoise = tansmitFunction(randomInputBits, signalVariance, SNR);
        bitErrorRateTrails(j) = calculateBitError(randomInputBits, outputWithNoise);
        
    end
    %Calculate the mean errorrate as the final result of certain SNR
    bitErrorRateOutput(i+1) = mean(bitErrorRateTrails); 
end


semilogy(SNRAxis, bitErrorRateOutput);
axis([0 50 -2 1])
title("Bit Error Rate vs SNR");
xlabel('E_{b}/N_{0}') ;
ylabel('P_{e}') ;


%---------------------------------------------------------------------
%Function Code
%---------------------------------------------------------------------



function outputWithNoise = tansmitFunction(input, inputVariance, SNR)
    % Generate equal number of noise samples
    % The generated noise should have normal distribution with zero mean and unit variance (use the function randn in MATLAB).
    mean = 0;
    numberOfNoiseBits = 1024;
    
    % Change the noise variance with respect to SNR (signal to noise ratio) value
    % Use SNR value to generate noise variance
    % SNR (in dB) = 10log10 (S/N) where S is the Signal power (or variance) and N is the Noise power (or variance)
    noiseVariance = inputVariance / 10^(SNR/10);
    noiseSTD = sqrt(noiseVariance);
    noiseBits = noiseSTD.*randn(numberOfNoiseBits, 1) + mean;
    
    %Find output with noise data
    outputWithNoise = input + noiseBits;
    %Fix the threshold value as 0 (the transmitted data is +1 and -1, and 0 is the midvalue)
    outputWithNoise(outputWithNoise > 0) = 1;
    outputWithNoise(outputWithNoise < 0) = -1;
    
    
    
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
