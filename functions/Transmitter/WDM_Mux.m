function [Sout,Pout] = WDM_Mux(S,P)
%% WDM COMPLEX SIGNALS MODULATOR
% Generates a time-domain multi wavelength signal from a series of
% pulse shaped sequences (multiple channels) performing upsampling when
% required
%
% INPUTS:
% P
%  .Sys.Nch                  Number of WDM channels
%  .Sys.Chsp                 Channel spacing
%  .Sim.Or                   Oversampling ratio
%
%% WDM parameters
% P.Sys.Chsp                 WDM channel spacing (assumed uniform)
% P.Sys.Npol                 Number of polarisations
% P.Sys.Nch                  Number of channels
%
%% Optional
% P.Sys.Pch                 1xNch vector of powers per channel [dBm] (if single value assumes same power over all channels)
%                            % indexed as -(N-1)/2,-(N-1)/2+1,...,(N-1)/2
%
% OUTPUTS
%
% Pout                     Output parameter structure
% Sout.Et                  Output WDM signal in time domain
%
% Author: Gabriele Liga, January 2019

Pout = P;
Sout = S;

%% Prepare output
Sout.Et = zeros(Pout.Sys.Npol,size(S.Et,2)*(Pout.Sim.Fs/S.Fs));
if ~isequal(Pout.Sim.Fs,S.Fs)
    Sout.Fs = Pout.Sim.Fs;
    Sout = MakeTimeFrequencyArray(Sout);
end

%% WDM Grid
if mod(P.Sys.Nch,2)
    channel_no = -(P.Sys.Nch-1)/2:(P.Sys.Nch-1)/2;      % normalized center frequencies
else
    channel_no = [(-P.Sys.Nch/2:-1)+1/2,(1:P.Sys.Nch/2)-1/2];
end

F_Grid = P.Sys.Chsp*channel_no;                       % ideal frequency grid in Hz
dF = Pout.Sim.Fs/length(Sout.FF);                       % simulation frequency resolution
Pout.Sys.Fchan = round(F_Grid/dF)*dF;                     % rounded frequency grid in Hz

%% WDM channel loop
plch = 10.^((P.Sys.Pch-30)/10);  % vector of powers per channel in linear units [W]

for pp = 1:Pout.Sys.Npol
    for nn = 1:(P.Sys.Nch)
        ChIdx = (nn-1)*Pout.Sys.Npol+pp;

        % Resampling
        if ~isequal(Pout.Sim.Fs,S.Fs)
            Et_idx = interpft(S.Et(ChIdx,:),size(Sout.Et,2));
        else
            Et_idx = S.Et(ChIdx,:);
        end

        % Normalise input signal power
        Poutw = mean(abs(Et_idx).^2);
        Et_idx = Et_idx/sqrt(Poutw);

        % Upconversion
        CompExp = exp(-2j*pi*Pout.Sys.Fchan(nn)*Sout.TT);
        if length(plch) == P.Sys.Nch
            Sout.Et(pp,:) = Sout.Et(pp,:) + sqrt(plch(nn)/Pout.Sys.Npol)*Et_idx.*CompExp;
        else
            Sout.Et(pp,:) = Sout.Et(pp,:) + sqrt(plch/Pout.Sys.Npol)*Et_idx.*CompExp;
        end
    end
end
