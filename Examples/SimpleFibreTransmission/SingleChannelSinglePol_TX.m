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
[Smux,P] = WDM_Mux(S,P);

%% Fibre propagation
fprintf('Fibre channel propagation...\n');
for ss = 1:P.Link.nspans
    fprintf('Span No %i\n',ss);
    [Smux,P] = FibreChannel(Smux,P);
    if P.Link.EDFA
        Smux = EDFA(Smux,P);
    end
end

%% Resampling and separate WDM channels
[Sdmux,P] = WDM_DeMux(Smux,P);

%% Matched filter
Smf = MatchedFilter(Sdmux,P);

%% Phase Recovery
Sout = ConstSync(Smf,P); 

%% Calculate metrics before demapping
R = CalcPerfMetric(R,Sout,P,'SNR');

heatscatterplot(Sout.RxSym.EstimatedSymbols);

