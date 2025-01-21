%% Run simple fiber transmission simulation
%--------------------------------------------------------------------------
% Mar. 2019, Sebastiaan Goossens
% Jan. 2021, Yunus Can Gültekin
%--------------------------------------------------------------------------
clear;
close all;

%% Parameters Initilization (S and P)
Parameters_SingleChannelSinglePol;
P = ParseParameters(P);

%% Convert optical units to SI units
P = Opt2SI(1,P);

%% Generate constellation and labeling
P = ModGen(P);

%% Generate information bits
[S, P] = InfoGen(P);

%% Generating random sequence of symbols from selected constellation
S = SeqGen(S,P);

%% Pulse shaping
S = PulseShape(S,P);

%% Resampling, combine WDM channels and scale power
[S,P] = WDM_Mux(S,P);

%% Fibre propagation
fprintf('Fibre channel propagation...\n');
for ss = 1:P.Link.nspans
    fprintf('Span No %i\n',ss);
    [S,P] = FibreChannel(S,P);
    if P.Link.EDFA
        S = EDFA(S,P);
    end
end

%% Resampling and separate WDM channels
[S,P] = WDM_DeMux(S,P);

%% Dispersion compensation
S = EDC(S,P);

%% Matched filter
S = MatchedFilter(S,P);

%% Phase Recovery
S = ConstSync(S,P); 

%% Calculate metrics before demapping
R = CalcPerfMetric(R,S,P,'SNR');

heatscatterplot(S.RxSym.EstimatedSymbols);

