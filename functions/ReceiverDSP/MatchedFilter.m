function Sout = MatchedFilter(S,P)
%% This function filters the received signal and downsamples it
%
% Input:
% S     -Signal structure
% P     -Filter structure
%
% Output:
% Sout.Et    -The downsampled version of the input signal
%
% Author: Sebastiaan Goossens, March 2019

%%
S.Fs=P.Sys.Rs*P.Rx.Ns;
S.Ns=P.Rx.Ns;
P.Filter.shape=P.Tx.Filter.shape;               % Pulse shaping. Takes as default 'RRC' Must specify ['RRC', 'Rect', 'Gauss','Generic']
P.Filter.type=P.Tx.Filter.type;                 % Must specify ['FIR'] takes as default 'Ideal'
P.Filter.BW=P.Sys.Rs;                          % Bandwidth (3dB cutoff freq. for Gauss) of the filter

% Select filter type
if isequal(P.Tx.Filter.shape, 'RRC')
    P.Filter.RRCrolloff=P.Tx.Filter.RRCrolloff; % Root-raise cosine roll-off
end
if isfield(P.Tx.Filter, 'taps')                   % Vector containing custom taps for FIR filter
    P.Filter.taps=P.Tx.Filter.taps;
end

Sout=Filter(S,P);
Sout.RxSym.EstimatedSymbols=downsample(Sout.Et.',P.Rx.Ns).';
end