function Pout = Initialize_FEC(P)
%--------------------------------------------------------------------------
% This function updates FEC-related parameters.
% Jan. 2021, Created: Yunus Can Gültekin
% Feb. 2022, Modified: LDPC Added. Kaiquan Wu
%--------------------------------------------------------------------------
% Input
% P: System Paramaters
% S: Signal Paramaters
%--------------------------------------------------------------------------
% Output
% P: System Paramaters for FEC struct
%    - .n, codeword length
%    - .k, the information bit length for each codeword
%    - .p, parity bit length for each codeword
%    - .Nblock, the number of FEC codewords
%    - .H, the parity check matrix (LDPC)
%    - ...
%--------------------------------------------------------------------------

Pout = P;
Pout.m = log2(Pout.M);

switch lower(P.FEC.CodeType)

    %======== LDPC FEC
    case {'dvbs2_ldpc','wifi_ldpc'}

        switch lower(P.FEC.CodeType)
            case 'dvbs2_ldpc' % DVB-S2 LDPC FEC
                Pout.FEC.n = 64800;
                Pout.FEC.H = dvbs2ldpc(Pout.FEC.R);
            case 'wifi_ldpc' % 802.11 LDPC FEC
                Pout.FEC.n = 648; % codeword length can be 648, 1296, 1944
                Pout.FEC.H = wlan11ParityCheckMatrices(Pout.FEC.R ,Pout.FEC.n);
        end

        [Pout.FEC.p,Pout.FEC.n] = size(Pout.FEC.H); % n-k parity bits, and n codeword length
        Pout.FEC.k = Pout.FEC.n - Pout.FEC.p;      % k is also equal to (R)*(n)
        Pout.FEC.LDPCenccfg = ldpcEncoderConfig(logical(sparse(Pout.FEC.H)));
        try
            Pout.FEC.LDPCdec = comm.gpu.LDPCDecoder(Pout.FEC.H);
        catch
            Pout.FEC.LDPCdeccfg = ldpcDecoderConfig(Pout.FEC.LDPCenccfg);
        end

        %======== PC-FEC
    case 'pc'
        Pout.FEC.v=8;
        Pout.FEC.t=3;
        Pout.FEC.s=0;
        Pout.FEC.itermax=10;
        Pout.FEC.nc=2^(Pout.FEC.v)-1;                                      % original component block length
        Pout.FEC.kc=2^(Pout.FEC.v)-Pout.FEC.v*Pout.FEC.t-1;                % original component information length
        Pout.FEC.ns=2^(Pout.FEC.v)-1-Pout.FEC.s;                           % Shorthened component block length
        Pout.FEC.ks=2^(Pout.FEC.v)-Pout.FEC.v*Pout.FEC.t-1-Pout.FEC.s;     % Shorthened component information length
        Pout.FEC.n = Pout.FEC.ns^2;                                          % product code block length
        Pout.FEC.k = Pout.FEC.ks^2;                                        % product code information length
end

if isfield(P.FEC,'Nfecframe') && isfield(P.Tx, 'Nsym')
    Pout.Tx.Nsym      = ceil(Pout.FEC.n * P.FEC.Nfecframe / Pout.m);
elseif isfield(P.FEC,'Nfecframe') && ~isfield(P.Tx, 'Nsym')
    Pout.Tx.Nsym      = ceil(Pout.FEC.n * P.FEC.Nfecframe / Pout.m);
elseif ~isfield(P.FEC,'Nfecframe') && isfield(P.Tx, 'Nsym')
    Pout.FEC.Nfecframe = ceil((P.Tx.Nsym*Pout.m)/Pout.FEC.n);            % the number of FEC codewords are obtained given the desired symbol numbers
else
end

if Pout.FEC.Nfecframe<1
    Pout.FEC.Nfecframe=1;
    warning('Instead of a # of FEC frames, a # of TX symbols that corresponds to <1 FEC frames is specifed. Now adjusted to TX 1 frame.')
end

Pout.FEC.Nbit  = Pout.FEC.Nfecframe * Pout.FEC.k;                          % number of information bits per channel per polarization
Pout.FEC.Ncbit = Pout.FEC.Nfecframe * Pout.FEC.n;                          % Number of coded bits per channel per polarization

end
