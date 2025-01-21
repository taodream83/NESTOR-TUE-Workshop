function Sout = SpectrumAnalyser(S,P)
%% Spectrum analyser function
% Plot and return abs squared value of the spectrum of the input signal (in time-domain) in log units
% If S.Et is a multidimensional signal the spectrum is the sum of the spectra over all dimensions
% Spectrum plot is shown in log units (dB)
%
% Inputs:
% S          - Signal structure
%      .Et       - Field
%      .Fs       - Signal sampling frequency
% P          - Filter structure
%   P.Filter.type -> Type of the filter: 'Ideal' or 'FIR'
%   (Ideal frequency implementation or FIR implementation with finite number of taps
%
%   P.Filter.shape -> Shape of the filter: 'RRC' (Root Raised CoSe),
%   'Rect' (Rectangular filter) or 'Gaussian' (Gaussian Filter)
%
%   P.Filter.FF -> Frequency Vector (not fftshifted, aligned with fft output)
%
%   P.Filter.BW -> Bandwidth of the Filter (for Gaussian, 3dB cutoff)
%
% Returns:
% Sout
%       .FF         - Frequency array [Hz]
%       .TT         - Time array [s]
%
% Author: Gabriele Liga, February 2019
Sout=MakeTimeFrequencyArray(S);
FF=fftshift(Sout.FF);
Nt=length(Sout.TT);
%P.Meter.Res
% Calculate normalised power spectral density
Sout.Et=fftshift(sum(abs(fft(S.Et,[],2)).^2,1));


if isfield(P, 'Meter') && isfield(P.Meter,'FreqRes')
    Res=P.Meter.FreqRes;
    dF=S.Fs/length(FF);
    NN= 2*round(Res/(2*dF));
    b= ones(1,NN);
    Sout.Et = filter(b, 1, [Sout.Et(NN:-1:1) Sout.Et Sout.Et(Nt:-1:Nt-NN-1)]);
    Sout.Et=Sout.Et(NN+NN/2+1:Nt+NN+NN/2);
end

if ~(isfield(P, 'Meter') && isfield(P.Meter,'FreqScale'))
    P.Meter.FreqScale='GHz';
end
switch P.Meter.FreqScale
    case 'Hz'
        flabel='Frequency [Hz]';
    case 'GHz'
        FF=FF*1e-9;
        flabel='Frequency [GHz]';
    case 'THz'
        FF=FF*1e-12;
        flabel='Frequency [THz]';
end

%N=size(S.Et,1);
Sout.Et=Sout.Et/max(Sout.Et);
figure
plot(FF,10*log10(Sout.Et),'k','LineWidth',1.2);

if isfield(P, 'Meter') && isfield(P.Meter,'Range')
    ylim([P.Meter.Range/10-P.Meter.Range, P.Meter.Range/10]);
    xlim([FF(1) FF(end)]);
    grid on;
    xlabel(flabel);
    ylabel('Normalised Power Spectral Density [dB/Hz]');
end