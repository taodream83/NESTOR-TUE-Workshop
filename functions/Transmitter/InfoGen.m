function [Sout,Pout] = InfoGen(P)
%--------------------------------------------------------------------------
% This function creates the information bits to be transmitted.
% Jan. 2021, Created: Yunus Can Gültekin
% Apr. 2021, Modified: Yunus Can Gültekin
% Feb. 2021, Modified: LDPC and PAS enabled. Kaiquan Wu
%--------------------------------------------------------------------------
% Input
% P: System Paramaters
% S: Signal Paramaters
%--------------------------------------------------------------------------
% Output
% P: System Paramaters
%    - .BitSeq/InfoBits, the information bits
% S: Signal Paramaters
%    - .Tx.Nsym/FEC.Ncsym, update the numner of coded syms per channel per polarization
%    - .FEC.Nblock, the number of FEC codewords
%    - .SHA, shaping related struct if do PAS
%--------------------------------------------------------------------------

Pout = P;
%---- Existing Parameters
Npol = P.Sys.Npol;
Nchn = P.Sys.Nch;
m_ND = P.m;
%----

% In the binary matrices outputted here, each row represents a polarization
% which consists of Nsym ND symbols, each addressed by m_ND bits

if strcmpi(P.FEC.CodeType,'none') && strcmpi(P.SHA.Type,'none')
    %% Unshaped & uncoded transmission, i.e., all m_ND bits of a ND-symbol are info
    Sout.BitSeq      = randi([0,1], 2*Nchn*Npol/P.N, P.Tx.Nsym*m_ND);

elseif ~strcmpi(P.FEC.CodeType,'none') && strcmpi(P.SHA.Type,'none')
    %% Unshaped & coded transmission, i.e., Rc*m_ND bits of a ND-symbol are info
    dimNum = P.Sys.Npol*P.Sys.Nch;     % Number of dimension strem
    Sout.InfoBits = randi([0,1],dimNum,P.FEC.Nbit); % Information bit sequence

elseif strcmpi(P.FEC.CodeType,'none') && ~strcmpi(P.SHA.Type,'none')
    %% Shaped & uncoded transmission, i.e., 2*(Rs+1) bits of a 2D-symbol are info (Rs: shaPg rate in bit/1D)

    if strcmpi(P.SHA.Type,'ideal')
        % Do nothing, symbols are drawn randomly from a predefined PMF. No
        % actual information bits are transmitted.
    else
        if mod(2*P.Tx.Nsym,P.SHA.N)
            warning('Adjusted P.Tx.Nsym to fit chosen blocklength');
        end
        Pout.Tx.Nsym = P.SHA.N * ceil(2*P.Tx.Nsym/P.SHA.N) / 2;

        if Pout.SHA.Kshp >= Pout.SHA.Ks
            Sout.InfoBits = randi([0,1], Npol*Nchn, Pout.Tx.Nsym*2*(Pout.SHA.Rs+1) );
        else
            error('Amplitude shaper configuration does not lead to a shaping rate that satisfies the target.');
        end
        Pout.SHA.Nblocks = 2*Pout.Tx.Nsym/Pout.SHA.N; % the total number of CCDM amplitude blocks
    end
elseif ~strcmpi(P.FEC.CodeType,'none') && ~strcmpi(P.SHA.Type,'none')
    %% Shaped & coded transmission, i.e., 2*(Rs+gamma) bits of a 2D-symbol are info (`gamma' requires some knowledge on PAS)

    % Bit string size for each FEC block
    Pout.SHA.NgammaBitFEC = 2*Pout.FEC.n/Pout.m-(Pout.FEC.n-Pout.FEC.k); % the number of info bits for sign bits per FEC frame
    Pout.SHA.NampsFEC = (Pout.FEC.n-Pout.FEC.k)+Pout.SHA.NgammaBitFEC; % the number of DM output amplitudes per FEC frame
    Pout.SHA.NampBitFEC = Pout.SHA.Rs*Pout.SHA.NampsFEC; % the number of info bits for DM input bits per FEC frame
    Pout.SHA.NinfoBitFEC = Pout.SHA.NgammaBitFEC+Pout.SHA.NampBitFEC; % the number of total info bits per FEC frame
    % the block number for FEC codewords and shaping blocks
    Pout.FEC.Nblock = 2*Pout.Tx.Nsym/Pout.SHA.NampsFEC; % 2 due to QAM
    if round(Pout.FEC.Nblock)==Pout.FEC.Nblock % checking if integer
    else
        error('Blocklength mismatch b/w SHA and FEC.')
    end
    Pout.SHA.Nblocks = Pout.SHA.NampsFEC*Pout.FEC.Nblock/Pout.SHA.N; % the total number of DM amplitude blocks

    % the shaping blocklength and FEC codeword length should generate a
    % integer number of amplitude blocks
    if round(Pout.SHA.Nblocks)~=Pout.SHA.Nblocks
        error('Shaping blocklength P.SHA.N is not compatiable with the FEC codeword length');
    end

    if mod(2*P.Tx.Nsym,P.SHA.N)
        warning('Adjusted P.Tx.Nsym to fit chosen blocklength');
    end
    Pout.Tx.Nsym = Pout.FEC.Nblock * Pout.SHA.NampsFEC/2;

    Sout.InfoBits = randi([0,1], Npol*Nchn, Pout.FEC.Nblock*Pout.SHA.NinfoBitFEC);


end

end

