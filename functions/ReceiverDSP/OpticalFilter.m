function [Sout] = OpticalFilter(S,P)
%% This function applies an optical filter to the received signal
%
% Input:
% S     -Signal structure
% P     -Filter structure
%
% Output:
% Sout.Et    -The optical filtered version of the input signal
%
% Author: Vinicius Oliari, February 2022
%%
S.Fs = P.Sim.Fs;
if isfield(P.Sim,'Ns')
    S.Ns = P.Sim.Ns;
else
    S.Ns = P.Rx.Ns;
end

if isfield(P.Sys,'OptFilter')
    if ~isfield(P.Sys.OptFilter,'shape')
        P.Filter.shape = 'Brickwall';
    else
        P.Filter.shape = P.Sys.OptFilter.shape;               % Pulse shaping. Takes as default 'RRC' Must specify ['RRC', 'Brickwall', 'Gauss', 'Generic']
        % Define roll-off in case of RRC
        if isequal(P.Sys.OptFilter.shape, 'RRC')
            P.Filter.RRCrolloff = P.Sys.OptFilter.RRCrolloff; % Root-raise cosine roll-off
        end
    end
else
    P.Filter.shape = 'Brickwall';
end

P.Filter.type  = 'Ideal';                 % Here, optical filters consider 'Ideal' filtering as default
P.Filter.BW    = P.Sys.Chsp;                          % Bandwidth (3dB cutoff freq. for Gauss) of the filter

Sout=Filter(S,P);
Sout.RxSym.EstimatedSymbols = downsample(Sout.Et.',P.Rx.Ns).';

end

