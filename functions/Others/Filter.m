function Sout = Filter(S,P)
% Filter
%
% This function filters an input signal according to a filter structure
% P
%
% %%%% Remarks %%%%
%
% % Necessary variables for all type of filters
%
%   P.Filter.type -> Type of the filter: 'Ideal' or 'FIR'
%   (Ideal frequency implementation or FIR implementation with finite number of taps
%
%   P.Filter.shape -> Shape of the filter: 'RRC' (Root Raised Cosine),
%   'Rect' (Rectangular filter) or 'Gaussian' (Gaussian Filter)
%
%   P.Filter.FF -> Frequency Vector (not fftshifted, aligned with fft output)
%
%   P.Filter.BW -> Bandwidth of the Filter (for Gaussian, 3dB cutoff)
%
% % Necessary variables for Ideal filtering
%
% % % RRC
%   P.Filter.RRCrolloff -> Roll-Off parameter of the Root Raised Cosine Filter
%
% % Necessary variables for FIR filtering
%
% % % RRC
%   P.Filter.RRCrolloff -> Roll-Off parameter of the Root Raised Cosine Filter
%   P.Filter.Ntaps -> Number of taps for the FIR filter
%
% % % Rect
%   P.Filter.Ntaps -> Number of taps for the FIR filter
%
% % % Generic
%   P.Filter.genFIR -> generic filter with an arbitrary number of taps
%
% % Frequency Normalzation
%
% As the frequencies are normalized in the beginning of the function, new
% variables added that are frequency related will also need to be
% normalized. In addition, many values and calculations slightly simplified
% due to the normalization.
%
% %%%% INPUTS %%%%%
%
% S: Input Signal Structure -> S.Et is the signal (NxT1) to be filtered,
% where T1 is the number of time samples
%
% P: Filter Structure -> Structure containing the parameters of the
% filter
%
% %%%% OUTPUT %%%%%
%
% Sout: Output Signal Structure -> Sout.Et is the desired filtered
% output signal (NxT1)
% Author: Vinicius Oliari, 29/01/2019

if P.Flags.DoublePrecision
    precision = @(x)double(x);
else
    precision = @(x)single(x);
end
if P.Flags.GPU
    dataFunc = @(x)gpuArray(precision(x));
else
    dataFunc = @(x)precision(x);
end

% Normalization of the frequencies (better accuracy for high frequencies)

Sout=MakeTimeFrequencyArray(S);
FF = Sout.FF; % Make a short variable name
maxF = max(abs(FF)); % Maximum frequency value
FF       = FF/maxF; % Normalize frequency vector
BW = P.Filter.BW /maxF; % Normalize filter bandwidth/cutoff frequency



switch P.Filter.type
    case 'Ideal'
        switch P.Filter.shape
            case 'RRC'

                FTSig = fft(S.Et,[],2); % Input signal in the Fourier domain

                rolloff = P.Filter.RRCrolloff; % Make a shorter variable name - RollOff of the RRCOS

                %%% Creates the RC filter
                SFreq = (1-rolloff)*BW/2; % Beggining of rolloff frequencies
                EFreq = (1+rolloff)*BW/2; % End of rolloff frequencies

                IndMid =  abs(FF)<= SFreq ; % Frequency positions where the gain is 1
                IndEdg = find( abs(FF)> SFreq & abs(FF)<=EFreq ); % Frequency positions of the rolloff

                RCFilt = zeros(1,length(FF)); % Initialize Filter vector
                RCFilt(IndMid) = 1; % Attribute 1 to the center frequencies (flat gain region)

                % Construct the rolloff region in the positions IndEdg
                RCFilt(IndEdg) = 0.5*(1+cos((pi/(rolloff*BW))...
                    *(abs(FF(IndEdg))-(1-rolloff)*BW/2)));
                %%%

                RRCFilt = sqrt(RCFilt); % Takes the squareroot of the RC filter (creating the RRC)

                FOutSig = FTSig.*RRCFilt; % Filtered signal (Fourier Domain)
                Sout.Et = ifft(FOutSig,[],2); % Comes back to the time domain

            case 'Brickwall'

                FTSig = fft(dataFunc(S.Et),[],2); % Input signal in the Fourier domain

                %%% Creates the Rect filter
                IndMid =  abs(FF)<= BW/2; % Frequency positions where the gain is 1

                %RectFilt = zeros(1,length(FF),'gpuArray'); % Initialize Filter vector
                RectFilt = zeros(1,length(FF)); % Initialize Filter vector

                RectFilt(IndMid) = 1; % Attribute 1 to the center frequencies (flat gain region)
                %%%

                FOutSig = FTSig.*RectFilt; % Filtered signal (Fourier Domain)
                Sout.Et = gather(ifft(FOutSig,[],2)); % Comes back to the time domain

            case 'Gauss'

                FTSig = fft(S.Et,[],2); % Input signal in the Fourier domain

                Freq3dB = BW; % Attributes the cutoff frequency (3dB) to P.Filter.BW

                % Computes the variace of the gaussian function with
                % respect to the frequency where the response drops 3dB
                sigmafsq = -(Freq3dB^2)/(2*log(0.5));

                GFilt = exp(-FF.^2/(2*sigmafsq)); % Creates the gaussian filter

                FOutSig = FTSig.*GFilt; % Filtered signal (Fourier Domain)
                Sout.Et = ifft(FOutSig,[],2); % Comes back to the time domain

            otherwise
                % Error Message - none of the available filter shapes
                disp('Error: Invalid Filter Shape. Try RRC or Rect.');
        end
    case 'FIR'
        switch P.Filter.shape
            case 'RRC'

                rolloff = P.Filter.RRCrolloff; % Make a shorter variable name - RollOff of the RRCOS
                Ntaps = P.Filter.Ntaps;  % Make a shorter variable name - Number of Taps
                % Converts the two-side bandwidth to one-side bandwidth
                btratio   = BW/2;

                FIR_N = GenRRCFIR(Ntaps,rolloff,btratio); % Generates the RRC pulse (time domain)

                Sout.Et = Filter_FIR(S.Et,FIR_N); % Applies the filter to the Signal (Circular Conv)

            case 'Brickwall'
                % As this function uses the degeneration of the RRC filter
                % to generate the retangular filter, we set rolloff=0
                rolloff = 0;
                Ntaps = P.Filter.Ntaps;  % Make a shorter variable name - Number of Taps
                % Converts the two-side bandwidth to one-side bandwidth
                btratio   =BW/2;

                FIR_N = GenRRCFIR(Ntaps,rolloff,btratio); % Generates the Sc pulse (time domain)

                Sout.Et = Filter_FIR(S.Et,FIR_N); % Applies the filter to the Signal (Circular Conv)

            case 'Gauss'

                Freq3dB = BW; % Attributes the cutoff frequency (3dB) to P.Filter.BW

                % Computes the variace of the gaussian function with
                % respect to the frequency where the response drops 3dB
                sigmafsq = -(Freq3dB^2)/(2*log(0.5));

                % Determine the gaussian variance in the time domain
                sigmatsq   = 1/(pi*pi*sigmafsq);

                % Sample time (normalized in the frequency normalization)
                Ts = 1;

                % Creates a time vector with the correct amount of taps
                if mod(P.Filter.Ntaps,2)==1
                    tv = (-floor(P.Filter.Ntaps/2):floor(P.Filter.Ntaps/2))*Ts;
                else
                    tv = (-floor(P.Filter.Ntaps/2):floor(P.Filter.Ntaps/2)-1)*Ts;
                end

                % Creates the gaussian vector
                FIR_N  = (Ts/sqrt(2*pi*sigmatsq))*exp(-(tv.^2)/(2*sigmatsq));

                Sout.Et = Filter_FIR(S.Et,FIR_N); % Applies the filter to the Signal (Circular Conv)

            case 'Generic'

                FIR_N = P.Filter.genFIR; % Make a shorter variable name - The generic FIR filter

                Sout.Et = Filter_FIR(S.Et,FIR_N); % Applies the filter to the Signal (Circular Conv)

        end
    otherwise
        % Error Message - none of the available filter shapes
        disp('Error: Invalid Filter Type. Try Ideal or FIR.');
end

end


function [OSig] = Filter_FIR(InSig,FIR_N)

[M,N] = size(InSig); % Obtain the signal dimensions

OSig = zeros(M,N); % Initialize the Output Signal

% Applies a loop in every row dimension - needed since the
% function cconv is only applied to vectors
for kk = 1:M
    % Makes a circular convolution of the signal with the
    % FIR filter
    OSig(kk,:) = cconv(InSig(kk,:),FIR_N,N);
end

end

function [RRC_FIR] = GenRRCFIR(Ntaps,rolloff,bwratio)

Ts = 1/(bwratio); % Normalized sample time

% Creates a time vector with the correct amount of taps
if mod(Ntaps,2)==1
    tv = (-floor(Ntaps/2):floor(Ntaps/2))*bwratio;
else
    tv = (-floor(Ntaps/2):floor(Ntaps/2)-1)*bwratio;
end

% Generates the RRC filter coefficients
RRC_FIR = (1/Ts)*(sin(pi*tv*(1-rolloff))+4*rolloff*tv.*cos(pi*tv*(1+rolloff)));
RRC_FIR = RRC_FIR./(pi*tv.*(1-(4*rolloff*tv).^2));

% Determine the points where the previous eq. was not defined
RRC_FIR(tv==0) = (1/Ts)*(1+rolloff*(4/pi-1));
RRC_FIR(abs(tv)==1/(4*rolloff)) = rolloff/(Ts*sqrt(2))*...
    ((1+2/pi)*sin(pi/(4*rolloff))+(1-2/pi)*cos(pi/(4*rolloff)));

end
