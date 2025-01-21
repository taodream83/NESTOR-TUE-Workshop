function [Sout,Pout] = SHA_DEC(S,P)
%--------------------------------------------------------------------------
% This function implements amplitude deshaping for PAS.
% Mar. 2022, Created. CCDM deshaping realized, ESS to be implemented. Kaiquan Wu
%--------------------------------------------------------------------------
% Input
% P: System Paramaters - No updates
%--------------------------------------------------------------------------
% S: Signal Paramaters
%    - .RxSym.EstimatedBitsHD, the recovered information bit sequences
%--------------------------------------------------------------------------

Sout = S;
Pout = P;

if ~strcmpi(P.SHA.Type,'none')
    %% Shaped transmission
    % get the 1D PAM amplitude set
    P_1D = P.P_1D;
    [~,idx_sort] = sort(P_1D.X);
    [~,ampIdxSet] = sort(idx_sort);

    % only decode dimensions from channel of interests
    pnt_dim = 1:P.Sys.Npol*length(P.Rx.CUT);
    DesInfBitSeq = zeros(length(pnt_dim),Pout.FEC.Nblock*Pout.SHA.NinfoBitFEC);
    for rowind = pnt_dim % each row is a polarization
        % extract the recovered info bits for sign bits
        FECrecBits = reshape(Sout.DecInfBitSeq(rowind,:),[],Pout.FEC.Nblock);
        recGammaBits = reshape(FECrecBits(1+Pout.SHA.NampsFEC*(P_1D.m-1):end,:),1,[]);
        % bits to amplitude
        recAmpBits = reshape(FECrecBits(1:Pout.SHA.NampsFEC*(P_1D.m-1),:),P_1D.m-1,[]);
        recAmpIdx = ampIdxSet(bi2de(recAmpBits.','left-msb')+sqrt(Pout.M)/2+1)-sqrt(Pout.M)/2-1;
        % Amplitude deshaping
        if strcmpi(P.SHA.Type,'ccdm')     % CCDM
            recAmpInfoBits = zeros(Pout.SHA.Nblocks*Pout.SHA.Kshp,1);
            for blk = 1:Pout.SHA.Nblocks
                pntBit = (blk-1)*Pout.SHA.Kshp+1:blk*Pout.SHA.Kshp;
                pntAmp = (blk-1)*Pout.SHA.N+1:blk*Pout.SHA.N;
                recAmpInfoBits(pntBit) = CCDM_decode(recAmpIdx(pntAmp).',Pout.SHA.N,Pout.SHA.Kshp,Pout.SHA.CC_Ni);
            end
        elseif strcmpi(P.SHA.Type,'ess') % ESS
            recAmpInfoBits = zeros(Pout.SHA.Nblocks*Pout.SHA.Kshp,1);
            for blk = 1:Pout.SHA.Nblocks
                pntBit = (blk-1)*Pout.SHA.Kshp+1:blk*Pout.SHA.Kshp;
                pntAmp = (blk-1)*Pout.SHA.N+1:blk*Pout.SHA.N;
                recAmpInfoBits(pntBit) = ESS_decode(P.SHA.Trellis,P.SHA.TrellisData,recAmpIdx(pntAmp).');
            end

        elseif strcmpi(P.SHA.Type,'bess') % B-ESS
            recAmpInfoBits = zeros(Pout.SHA.Nblocks*Pout.SHA.Kshp,1);
            for blk = 1:Pout.SHA.Nblocks
                pntBit = (blk-1)*Pout.SHA.Kshp+1:blk*Pout.SHA.Kshp;
                pntAmp = (blk-1)*Pout.SHA.N+1:blk*Pout.SHA.N;
                % recAmpInfoBits(pntBit) = ESS_decode(P.SHA.Trellis,P.SHA.TrellisData,recAmpIdx(pntAmp).');
                recAmpInfoBits(pntBit) = ESS_BP_Deshape(2*recAmpIdx(pntAmp)+1,P.SHA.Trellis,P.SHA.N,P.SHA.L,P.SHA.na,P.SHA.Lmant,P.SHA.Kshp);
            end

        else
            error('Invalid PS type specification!');
        end

        DesInfBitSeq(rowind,:) = [recAmpInfoBits(:).',recGammaBits]; % Each row represents a polarization

    end

    % For PAS, overwrite the recovered bits after FEC decoding
    if ~strcmpi(P.FEC.CodeType,'none')
        Sout.RxSym.EstimatedBitsHD = [];
        Sout.RxSym.EstimatedBitsHD = DesInfBitSeq;
    end
end


end
