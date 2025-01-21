function [Sout,Pout] = Resampling(S,P)

% Resampling - Vinicius Oliari - 24/01/2019
%
% This function resamples a signal S.Et that has P.Filter.Fsin sample
% rate, creating a new vector Sout.Et with P.Filter.Fsout
%
% %%%% INPUTS %%%%%
%
% S: Input Signal Structure -> Required Fields:
%
%      S.Et -> Input signal (Ndim x Nt matrix)
%
% P: P.Filter.Fsout -> Desired sampling rate
%      Pin.Filter.Fsin  -> Actual sampling rate
%
% %%%% OUTPUTS %%%%
%
% Sout: Output Signal Structure
%
% Pout: Output Parameters Structure
%

Sout=S;
Pout=P;
%Et=gpuArray(S.Et);
Sout.Et = interpft(S.Et.', (P.Filter.Fsout/P.Filter.Fsin)*length(S.Et)).';
% with the ratio P.Filter.Fsout/P.Filter.Fsin
Pout.Filter.Fsin=P.Filter.Fsout;      % Sets new input sampling frequency
Sout.Fs=Pout.Filter.Fsin;   % Sets current signal sampling frequency

end

