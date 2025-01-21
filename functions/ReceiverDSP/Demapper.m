function [Sout,Pout] = Demapper(S,P)
%% -------------------Demapper --------------------%%
% 1) MED symbol-by-symbol detection
% 2) HD demapping
% 3) SD demapping
% 4) LLR
% This function demap received symbol sequences only from the channel of interests.
% Jan. 2019, Created:  Bin Chen
% Apr. 2021, Modified: Yunus Can Gültekin
% Feb. 2022, Modified: Demapping for PAS enabled. Kaiquan Wu.
%--------------------------------------------------------------------------
% Input
% P: System Paramaters
%   - .Nsym, the number of symbols in one polarizations
%   - .Nbit, the number of transmitted bits over 1 tributary bit channel
%   - .Sys.Npol, the number of polarization
%   - .Sys.Nch, the number of Channel
%--------------------------------------------------------------------------
% S: Signal Paramaters
%    - .SymSeq, the tramsitted symbol sequences generated from input bits
%--------------------------------------------------------------------------

if ~strcmpi(P.FEC.CodeType,'none') && ~strcmpi(P.SHA.Type,'none')  % PAS case
    Sout = PASDemapper(S,P);
    Pout = P;
    return
end

Nsym=P.Tx.Nsym;

CUT=P.Rx.CUT;
Sout=S;
Pout=P;
Pout.PerfMetric.SER=zeros(1,length(CUT));
Pout.PerfMetric.BER=zeros(1,length(CUT));

%% Computes TX channel indices
if mod(P.Sys.Nch,2)
    CUTidx=P.Rx.CUT.'+(P.Sys.Nch-1)/2;
else
    CUTidx1=P.Rx.CUT(P.Rx.CUT<0).'+P.Sys.Nch/2;
    CUTidx2=P.Rx.CUT(P.Rx.CUT>0).'+P.Sys.Nch/2-1;
    CUTidx=[CUTidx1 ; CUTidx2];
end

TxChIdx = (repmat(1:P.Sys.Npol,length(P.Rx.CUT),1)+CUTidx*P.Sys.Npol).'; % Transmitted channel index
TxChIdx=TxChIdx(:);

%% MED symbol-by-symbol detection
% Initialise
Sout.RxSym.EstimatedIndicesHD=S.RxSym.EstimatedSymbols;
Sout.RxSym.EstimatedSymbolsHD=S.RxSym.EstimatedSymbols;

if P.N == 2   % Hard-decision on symbols
    LLR=zeros(length(TxChIdx),Nsym*P.m);
    for cc=1:length(TxChIdx)

        BB=repmat(P.C,Nsym,1);
        AA=repmat(S.RxSym.EstimatedSymbols(cc,:).',1,P.M);
        [~,EE]=min(abs(AA-BB),[],2);
        EE=EE.';               % In case EE is a (redundant) 3D array (dual-pol case)
        Sout.RxSym.EstimatedIndicesHD(cc,:)=EE;
        Sout.RxSym.EstimatedSymbolsHD(cc,:)=P.C(EE);
        %% Symbol-to-bit demapping
        EstBitsHD(cc,:)=reshape(de2bi(Sout.RxSym.EstimatedIndicesHD(cc,:).'-1,P.m,'left-msb').',Nsym*P.m,1).';

        %% Memoryless LLR calculation with signal-independent noise estimates
        llr = SoftDemapper(P.C.',S.SymSeq(TxChIdx(cc),:).',Sout.RxSym.EstimatedSymbols(cc,:).',P.Ik1,P.Ik0,P.PrX,P.N,0);
        LLR(cc,:) = reshape(llr,1,[]);
    end

elseif P.N == 4
    LLR=zeros(length(CUT),Nsym*P.m);
    for cc=1:length(CUT)
        %% Memoryless LLR calculation with signal-independent noise estimates
        llr = SoftDemapper(P.C.',S.SymSeq(TxChIdx(2*cc-1:2*cc),:).',Sout.RxSym.EstimatedSymbols(2*cc-1:2*cc,:).',P.Ik1,P.Ik0,P.PrX,P.N,0);
        LLR(cc,:) = reshape(llr,1,[]);
    end
    EstBitsHD = []; % Not supported yet
end


%% Assigns to output signal structure with or without depadding last bits for codeword bit matching
if  isfield(P,'FEC') && ~strcmpi(P.FEC.CodeType,'none')   % FEC case
    pad=Pout.FEC.pad;
    Sout.RxSym.EstimatedBitsHD=EstBitsHD(:,1:end-pad);
    Sout.RxSym.LLR= LLR(:,1:end-pad);
else
    Sout.RxSym.EstimatedBitsHD=EstBitsHD;
    Sout.RxSym.LLR= LLR;
end

end


%% Helper functions
function Sout = PASDemapper(S,P)
%--------------------------------------------------------------------------
% This function demapp the received symbol sequences into LLR
% Feb. 2022, Created: Kaiquan Wu
%--------------------------------------------------------------------------
% Input
% P: System Paramaters - No updates
%--------------------------------------------------------------------------
% S: Signal Paramaters
%    - .RxSym.LLR, the LLR sequences of received FEC codewords
%--------------------------------------------------------------------------
Sout = S;
% use the 1D PAM for demapping on each real dimension
% Extrac the 1D PAM signal struct
P_1D = P.P_1D;
P_1D.SHA.Type='none';P_1D.FEC.CodeType='none';
P_1D.PrX=[P_1D.SHA.PrA_1D(end:-1:1),P_1D.SHA.PrA_1D]./2;P_1D.PrX=P_1D.PrX.';P_1D.C = P_1D.X;

% Demapping in inphase/real of the received symbols
S_1D.RxSym.EstimatedSymbols = real(S.RxSym.EstimatedSymbols);
S_1D.SymSeq = real(S.SymSeq);
[S_1D,~] = Demapper(S_1D,P_1D);
realLLR = S_1D.RxSym.LLR;

% Demapping in imaginary/quadrature of the received symbols
S_1D.RxSym.EstimatedSymbols = imag(S.RxSym.EstimatedSymbols);
S_1D.SymSeq = imag(S.SymSeq);
[S_1D,~] = Demapper(S_1D,P_1D);
imagLLR = S_1D.RxSym.LLR;

% only demapp dimensions from channel of interests
pnt_dim = 1:P.Sys.Npol*length(P.Rx.CUT);
for rowind = pnt_dim
    % Extract the LLR from the receuved symbols on each dimension
    realLLRDim = reshape(realLLR(rowind,:),P_1D.m,[]);
    imagLLRDim = reshape(imagLLR(rowind,:),P_1D.m,[]);
    LLR= reshape([realLLRDim;imagLLRDim],P_1D.m,[]);
    % Categrize the LLRs
    signBitLLR = reshape(LLR(1,:),[],P.FEC.Nblock);
    parityBitLLR = signBitLLR(1:P.FEC.p,:);
    shapInfoBitLLR = signBitLLR(1+P.FEC.p:end,:);
    ampBitLLR = reshape(LLR(2:P_1D.m,:),[],P.FEC.Nblock);
    % Combined the LLRs
    Sout.RxSym.LLR(rowind,:) = reshape([ampBitLLR;shapInfoBitLLR;parityBitLLR],1,[]);
end
end

function LLR = SoftDemapper(C,X,Y,Ik1,Ik0,P,N,flag_demapper)
% This function computes L-values for a given constellation C based on the
% channel observations Y. The function returns a matrix of m times Ns LLRs.
% The input vectors are:
%
% C     : Column vector of M = 2^m constellation symbols
% X     : Column vector of complex transmitted symbols
% Y     : Column vector of complex received symbols
% Ik1   : M/2 times m matrix with pointers to the symbols in C labeled with
% a one in a given bit position.
% Ik0   : M/2 times m matrix with pointers to the symbols in C labeled with
% a zero in a given bit position.
% P     : Column vector of probabilities corresponding to C
% If flag_demapper == 0, then exact LLR calculation (Sum-Exp) is performed.
% If flag_demapper == 1, then approximated LLR calculation (Max-Log) is
% performed
%
% Alex Alvarado, November 2015
% Sebastiaan Goossens, 2020

Ns = size(Y,1);          % Number of received symbols
M = size(C,1);            % Constellation size
m = log2(M);             % Number of bits per symbol
VarZ = max(mean(var(Y-X)),1e-10);  % Noise variance
LLR = zeros(m,Ns);       % LLR output matrix

for k = 1:m
    dist1 = 0;
    dist0 = 0;
    for j = 1:N/2
        Ymat = repmat(Y(:,j),1,M/2); % For LLR calculation

        Cmat1 = repmat(C(Ik1(:,k),j).',Ns,1);
        Cmat0 = repmat(C(Ik0(:,k),j).',Ns,1);

        dist1 = abs(Ymat-Cmat1).^2 + dist1;
        dist0 = abs(Ymat-Cmat0).^2 + dist0;
    end

    Pmat1 = repmat(P(Ik1(:,k)).',Ns,1);
    Pmat0 = repmat(P(Ik0(:,k)).',Ns,1);

    if flag_demapper == 0 % Sum-Exp
        num = sum(Pmat1.*exp(-dist1/VarZ),2);
        den = sum(Pmat0.*exp(-dist0/VarZ),2);
    elseif flag_demapper == 1 %Max-log
        num = max(Pmat1.*exp(-dist1/VarZ),[],2);
        den = max(Pmat0.*exp(-dist0/VarZ),[],2);
    end
    LLR(k,:) = log(num) - log(den);
end
end