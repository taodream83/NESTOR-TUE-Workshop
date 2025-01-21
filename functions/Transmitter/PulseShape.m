function Sout = PulseShape(S,P)
%% This function produces a pulse shaped signal from an input sequence of symbols
%
% Input:
% S     -Signal structure
% P     -Filter structure
%
% Output:
% Sout    -Output signal structure
%
% Author: Gabriele Liga, Jan 2019


%% Upsampling
S.Et=upsample(S.SymSeq.',P.Tx.Filter.Ns).';
S.Fs=P.Sys.Rs*P.Tx.Filter.Ns;
% Passing parameters to Filter function
P.Filter.shape=P.Tx.Filter.shape;                        %Pulse shaping. Takes as default 'RRC' Must specify ['RRC', 'Rect', 'Gauss','Generic']
P.Filter.type=P.Tx.Filter.type;                      %Must specify ['FIR'] takes as default 'Ideal'
P.Filter.BW= P.Sys.Rs;                     %Bandwidth (3dB cutoff freq. for Gauss) of the filter

% Select filter type
if isequal(P.Tx.Filter.shape, 'RRC')
    P.Filter.RRCrolloff=P.Tx.Filter.RRCrolloff;                        % Root-raise cosine roll-off
end
if isfield(P.Tx.Filter, 'taps')                                 % Vector containing custom taps for FIR filter
    P.Filter.taps=P.Tx.Filter.taps;
end
Sout=Filter(S,P);

end