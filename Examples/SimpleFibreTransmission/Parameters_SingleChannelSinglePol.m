%% ------------------------ Initialization test --------------------- %%
% This file contains the parameters required for the corresponding example
% of a simulation of an optical fibre transmission system using the ICTlab optcomms toolbox


%% Constants and Flags                              %Include all the constants and flags to be used along the code
R = struct;                                        % Initialize results structure
P.Flags.GPU = 1;                                   % Enable GPU computing for SSFM
P.Flags.DoublePrecision = 1;                       % Define single or double precision datatype

%% Fiber                                            %Physical properties of the fiber
P.Fibre.PropModel='Manakov';                        % Fibre propagation model. Options are 'Manakov' for Manakov equation and 'NLSE' for the nonlinear Schroedinger equation.
P.Fibre.alpha       = 0.2;                          %Loss parameter  [dB/km]
P.Fibre.gamma       = 1.2;                          %Nonlinear parameter  [W/km]
P.Fibre.D           = 17;                           %Dispersion parameter at reference wavelength [ps/nm/km]


%% Link                                             %Physical properties of the links
P.Link.nspans       = 5;                            %Number of spans
P.Link.spanlength   = 80;                           %Span length [km]
P.Link.totlength  = P.Link.nspans*P.Link.spanlength;%Total length of the link [km]

% Optical amplifier                                 %Type of optical amplifier
P.Link.EDFA         = 1;                            %Erbium-Doped Fiber Amplifier flag
P.Link.Raman        = 0;                            %Raman Amplification flag
P.Link.NF           = 5;                            %Noise Figure [dB]


%% System                                           %System-transmition parameters
P.Sys.Nch        = 1;                            %Number of channels
P.Sys.Npol       = 1;                            %Number of polarization
P.Sys.Rs         = 60e9;                         %Symbol rate [sym/s]
P.Sys.Chsp       = 100e9;                         %Channel spacing [Hz]
P.Sys.Pch       =0;                              % 1xNch vector of powers per channel. If scalar all channels are assumed to have equal power [dBm]
P.Sys.lambda     =1550;                          % Reference lambda [nm]


%% Simulation                                       %Parameters of the simulation
P.Sim.Or=2;                                        % Oversampling ratio (with respect to simulation Nyquist BW)
P.Sim.Ns            =P.Sim.Or*ceil(P.Sys.Chsp/P.Sys.Rs)*P.Sys.Nch;            %Sampling rate (samples per symbol) [Sa/sym]
P.Sim.Fs            =P.Sim.Ns*P.Sys.Rs;              %Sampling rate [GSa/s]
P.Sim.dt            = 1/P.Sim.Fs;                   %Time step [ns]
P.Sim.SSFM.Type     = 'uniform';                    %Split-Step Fourier Method flavours. Options are ['uniform', 'log', 'adaptive']
P.Sim.SSFM.dz       =0.1;                           % Step size for 'uniform' implementation of SSFM [km]


%% Meter parameters (spectrum analyser)
P.Meter.FreqRes=1e9;            % Indicates frequency resolution of spectrum analyser [Hz]
P.Meter.Range=50;               % Dynamic range for spectrum analyser [dB]


%% Tx                                               %TX parameters
P.Tx.Nsym=1e4;                                     % Number of transmitted symbols over 1 polarisation
P.Tx.Filter.shape  = 'RRC';                        %Pulse shaping. Takes as default 'RRC' Must specify ['RRC', 'Brickwall', 'Gauss','Generic']
P.Tx.Filter.type   = 'Ideal';                      %Must specify ['FIR'] takes as default 'Ideal'
P.Tx.Filter.BW     = P.Sys.Rs;                     %Bandwidth (3dB cutoff freq. for Gauss) of the filter
P.Tx.Filter.RRCrolloff=0.01;                        % Root-raise cosine roll-off
P.Tx.Filter.Ns=2;                                   % Pulse shaping sampling frequency (samples/symbol)


%% Rx                                               %RX Parameters
P.Rx.Ns             = 2;                            % Samples/symbol
P.Rx.EDC.Custom     = 0;
P.Rx.ConstSyncMethod = 'xcorr';                       % Methods: 'xcorr'


%% FEC
% This substructure is for future developement
P.FEC.CodeType = 'none'; % 'none': no FEC, 'dvbs2_LDPC': DVB-S2 LDPC FEC, 'wifi_LDPC': IEEE 802.11 LDPC FEC, 'PC': Product codes


%% Signal structure %%%%%%                    %S is for signal (But this must be created in the transmission!)

%% Constellations Parameters
P.N           = 2;                            % Number of TX constellation real dimensions
P.M           = 16;                           % Cardinality in P.N (real) dimensions
P.ModFormat   = 'QAM';                        % Modulation format

