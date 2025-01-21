function Sout = MakeTimeFrequencyArray(S)
% Creates time and frequency vectors for current signal (S)
% Inputs:
% S          - Signal structure
%      .Et       - Field
%      .Fs       - Signal sampling frequency
%
% Returns:
% Sout
%       .FF         - Frequency array [Hz]
%       .TT         - Time array [s]
%
% Author: Gabriele Liga Jan 2019
Sout=S;
dT = 1/S.Fs;                                   % Temporal resolution (s)

%dF = 1/P.Sim.T;                                       % Frequency resolution (Hz)
Nt=size(S.Et,2);
T=dT*Nt;                                                 % Observation time
dF = 1/T;                                                  % Frequency resolution
Sout.TT = (0:Nt-1) * dT;                                 % Time array (s)
Sout.FF = [0:floor(Nt/2)-1,floor(-Nt/2):-1] * dF;        % Frequency array (Hz)

end

