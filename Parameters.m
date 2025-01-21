%% ------------------------ Parameters file --------------------- %%

% This file contains the comprehensive list of parameters required (and optional)
% for the simulation of an optical fibre transmission system using the ICTlab OpticalComm toolbox

%% WARNING: only some parameters are uncommented to provide a simple parameters' template for the unexperienced user.
%% Other parameters can be uncommented based on the specific simulation requirements (see toolbox documentation).

%% Constants and Flags                              %Include all the constants and flags to be used along the code
R = struct; % Initialize empty results structure
%P.Constants        = 0;
%P.Flags            = 1;
P.Flags.GPU = 1;      % Enable GPU computing for SSFM
P.Flags.DoublePrecision = 1;                       % Define single or double precision datatype

%% Fibre                                            %Physical properties of the fiber
P.Fibre.PropModel='Manakov';                          %Defines the propagation model adopted for the simulation of fibre propagation. Options are 'Manakov' for Manakov equation and 'NLSE' for the nonlinear Schroedinger equation. Default option is 'Manakov'.
P.Fibre.alpha       = 0.2;                           %Loss parameter  [dB/km]
%P.Fibre.beta2       = -21.67;                       %Second order dispersion coefficient [ps^2/km]
%P.Fibre.beta3       = 0;                            %Third order dispersion coefficient  [ps^3/km]
P.Fibre.gamma       = 1.2;                           %Nonlinear parameter  [W/km]
%P.Fibre.Biref       =1;                             % Flag determining whether or not simulate fibre birefringence. It applies a random rotation between the SOP of the optical field at the input and output of the fibre
%P.Fibre.PMD         = 0.1;                          % Polarization mode dispersion parameter [ps/(km)^1/2]
%P.Fibre.Lcorr       = 0.1;                          % Correlation length for SOP evolution (to be used in combination with PMD) [km]
P.Fibre.D           = 17;                            %Dispersion parameter at reference wavelength [ps/nm/km]
%P.Fibre.S           = 0;                            %Dispersion slope [ps/nm^2 km]
%P.Fibre.Cr         =0.008;                          % Raman gain slope [1/W/THz/km]
%P.Fibre.ISRS       =0;
%P.Fibre.Angles                                     %Stokes angles describing SOP evolution
%P.Fibre.DGD                                        %Vector of DGDs describing PMD evolution


%% Link                                             %Physical properties of the links
P.Link.nspans       = 1;                            %Number of spans
P.Link.spanlength   = 80;                           %Span length [km]
P.Link.totlength  = P.Link.nspans*P.Link.spanlength;%Total length of the link [km]
% Optical amplifier                                 %Type of optical amplifier
P.Link.EDFA         = 1;                            %Erbium-Doped Fiber Amplifier flag
%P.Link.Raman        = 0;                            %Raman Amplification flag
% P.Link.G            = 10;                           %Amplifiers' Gain [dB]
P.Link.NF           = 5;                            %Noise Figure [dB]

%% Laser                                            %Physical properties of the laser
%P.Laser.lambda      = 1550;                          %Laser wavelength [nm]

%% System                                           %System-transmition parameters
P.Sys.Nch        = 3;                            %Number of channels
P.Sys.Npol       = 2;                            %Number of polarization
P.Sys.Rs         = 32e9;                           %Symbol rate [Bd]
P.Sys.Chsp       = 50e9;                           %Channel spacing [Hz]
P.Sys.Pch       =20;                           % 1xNch vector of powers per channel [dBm]
P.Sys.lambda     = 1550;                        % Reference lambda [nm]

%% Optical filter parameters
P.Sys.OptFilter.shape   =  'Gauss' ;                       % Options for digital/optical filters
%P.Sys.OptFilter.BW;

%% Simulation                                       %Parameters of the simulation
P.Sim.Or=2;                                         % Oversampling ratio (with respect to simulation Nyquist BW)
P.Sim.Ns            =P.Sim.Or*ceil(P.Sys.Chsp/P.Sys.Rs)*P.Sys.Nch;            %Sampling rate (samples per symbol) [Sa/sym]
P.Sim.Fs            =P.Sim.Ns*P.Sys.Rs;              %Sampling rate [GSa/s]
%P.Sim.TT                                           %Time vector initialised by MakeTimeFrequenceArray.m [s]
%P.Sim.FF                                           %Frequency initialised by MakeTimeFrequenceArray.m [s]
P.Sim.dt            = 1/P.Sim.Fs;                   %Time step [ns]
%P.Sim.T                                            %Simulation time [s] (observation window)
P.Sim.SSFM.Type     = 'uniform';                    %Split-Step Fourier Method flavours. Options are ['uniform', 'log', 'adaptive']
P.Sim.SSFM.nsteps   =100;                           % Set number of steps per span instead of step size (for 'uniform', 'log' SSFM)
P.Sim.SSFM.dz       =0.05;
P.Sim.SSFM.PhiMax   =0.03;                          % Maximum nonlinear phase rotation allowed per step (only for 'adaptive' SSFM) [rad]

%% Optical filter parameters
P.Sys.OptFilter.shape = 'Brickwall';               % Pulse shaping. Takes as default 'RRC' Must specify ['RRC', 'Brickwall', 'Gauss','Generic']


%% Digital filter parameters
%P.Filter.Fsin;
%P.Filter.Fsout;
%P.Filter.shape  = 'RRC';
%P.Filter.BW;
%P.Filter.RRCrolloff=0.1;

%% Meter parameters (spectrum analyser)
%P.Meter.FreqScale='THz';         % Indicates unit for frequency axis in spectrum analyser.Options are 'Hz', 'GHz', or 'THz'
P.Meter.FreqRes=50e9;            % Indicates frequency resolution of spectrum analyser [Hz]
P.Meter.Range=50;               % Dynamic range for spectrum analyser [dB]


%% Tx                                               %TX parameters
P.Tx.Nsym=1e5;                                     % Number of transmitted symbols over 1 polarisation
%P.Tx.Nbit                                         % Number of transmitted bits over 1 tributary bit channel
P.Tx.Filter.shape  = 'RRC';                        %Pulse shaping. Takes as default 'RRC' Must specify ['RRC', 'Rect', 'Gauss','Generic']
P.Tx.Filter.type   = 'Ideal';                      %Must specify ['FIR'] takes as default 'Ideal'
P.Tx.Filter.BW     = P.Sys.Rs;                     %Bandwidth (3dB cutoff freq. for Gauss) of the filter
P.Tx.Filter.RRCrolloff=0.1;                        % Root-raise cosine roll-off
%P.Tx.Filter.taps;                                 % Vector containing custom taps for FIR filter
P.Tx.Filter.Ns=2;                                   % Pulse shaping sampling frequency (samples/symbol)

%% Rx                                               %RX Parameters
P.Rx.SaSym          = 2;                            % Samples/symbol
P.Rx.EDC.Custom     = 0;
P.Rx.CUT            = 1;                              % Index of channel under test, ranging from (1:P.Sys.Nch). 1 corespond to leftmost channel. Can be a vector if there are multiple interested WDM channels.
P.Rx.ConstSyncMethod = 'xcorr';                       % Methods: 'cxcov', 'rotnorm_average', 'rotnorm_cloud'
P.Rx.SymDisc=1e3;                                     % Number of discarded symbols before performance metric estimation to avoid cyclic effects.

% P.Rx.EDC.beta2      = P.Fibre.beta2;                % beta2 for EDC
% P.Rx.EDC.L          = P.Link.totlength;             % Length of compensation

% DBP parameters
P.DBP.Mode='Ideal';     % Options: 'Ideal'-ideally inverts channel,'Pragmatic'-inverts channel with specified complexity
%Pin.DBP.PMD=0;            % Flag to enable PMD inversion
%P.DBP.Type='uniform';         %  'log', 'adaptive'    Specifies DBP step distribution for pragmatic DBP mode
%P.DBP.PMD=0;                  % Flag to enable PMD inversion
%P.DBP.Nsteps=100;             % DBP steps/span for pragmatic DBP mode

%% FEC
% This substructure is for future developement
P.FEC.CodeType=3;                                    %0: no FEC 1:LDPC 2: HD-FEC 3: PC
P.FEC.v=8;
P.FEC.t=3;
P.FEC.s=0;
P.FEC.itermax=10;
%P.FEC.Nblock                  % Number of TX FEC blocks
%P.FEC.state                   % FEC state
%P.FEC.Encobj                  % Encoder object
%P.FEC.Decobj                  % Decoder object
%P.FEC.Ncbit        	       % Number of coded bits per channel per polarization
%P.FEC.Ncsym        	       % Number of coded syms per channel per polarization
%P.FEC.BitErr                  % Total number of bit errors in decoded sequence


%% Performance Metrics
%P.PerfMetric.BER;                 %  Uncoded BER
%P.PerfMetric.BERpost;
%P.PerfMetric.BERpre;
%P.PerfMetric.SER;
%P.PerfMetric.FER;
%P.PerfMetric.SNR;                 % Effective SNR [dB]
%P.PerfMetric.MI;                  % Mutual information [bits/2D]
%P.PerfMetric.GMI;
%P.PerfMetric.Q;                   % Q-factor [dB]

%% Signal structure                                   %S is for signal (But this must be created in the transmission!)
%S.BitSeq;                                              % Matrix of transmitted bit sequences per channel (S.BitSeqNt by P.Sys.Npol )
%S.SymSeq;                                              % Matrix of transmitted symbol sequences per channel
%S.Et;                                                  % Matrix of transmitted signal in time domain
%S.PSD;                                                 % Vector representing the discrete-time (transmitted) optical field power spectral density

% Constellation Parameters
%S.Fs;                                        Current signal sampling frequency
%S.Ns;                                        Current signal sampling frequency (samples/symbol)
%S.TT                                         Current signal time vector
%S.FF                                         Current signal frequency vector
S.N           = 2;                            % Number of TX constellation real dimensions
S.M           = 64;                           % Cardinality in S.N (real) dimensions
S.ModFormat   = 'QAM';                        % Modulation format
%S.C                                          % (S.N/2)xS.M vector of complex symbols representing the defined 2D or 4D constellation
%S.X                                          % Constellation
%S.m           = log2(S.M);                   % Bits/2Dsym


% S.Et                                                % Optical signal in the continuous-time domain
% S.SymSeq                                            % (P.Sys.NpolxP.Sys.Nch)xP.Tx.Nsym array of TX complex symbols
% S.TxIdx                                             % (P.Sys.NpolxP.Sys.Nch)xP.Tx.Nsym array of transmitted decimal indeces for TX constellation S.X
% S.BaseBand                                           % Indicates whether or not Et contains WDM channels at baseband
% S.DecInfBitSeq                                       % Decoded information bits


