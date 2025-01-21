function Sout = EDFA(S,P)

% Black box EDFA models gain, gain saturation and ASE
% The noise is split equally to both polarisation states
% Etout=Amplifier(dT,Etin,GdB,NFdB,PsatdBm)
% Inputs:
% S: input signal structure
%     -.Et   optical field


% P: input parameter structure
%    .Sys.Npol     - number of signal polarisations
%    .Sys.lambda;  - reference lambda              [m]
%    .Sim.Fs      - simulation sampling frequency [Hz]
%    .Link.G (optional)      - amplifier gain [linear units]
%    .Link.NF (optional)      - amplifier noise figure [linear units]

% Returns:
% SignalOut - output signal structure
% Author: Gabriele Liga, Jan 2019.

Nt=size(S.Et,2);
Np=P.Sys.Npol;
Sout=S;
fo=P.constant.c/P.Sys.lambda;

if isfield(P.Link, 'G')
    G=P.Link.G;     % Gain in linear units
else % If gain is not defined, EDFA will compensate for
    G=exp(P.Fibre.alpha*P.Link.spanlength);      % Gain in linear units
end

if isfield(P.Link, 'NF')
    NF=P.Link.NF;
    Nsp=(NF*G)/(2*(G-1));               % Calculate spontaneous emission factor from gain and noise figure
    Pase=2*Nsp*(G-1)*P.constant.h*fo*P.Sim.Fs;     % ASE power over simulation bandwidth over 2 polarisations (W)

    % White Gaussian noise zero mean and standard deviation equal to ASE
    % power split equally across the dimensions
    noiset = sqrt(Pase/4)*(randn(Np,Nt)+1i*randn(Np,Nt));
    Sout.Et = sqrt(G)*S.Et+noiset;
else
    Sout.Et = sqrt(G)*S.Et;
end

