function [Sout,Pout] = SHA_ENC(S,P)
%--------------------------------------------------------------------------
% This function implements amplitude shaping.
% Jan. 2021, Created: Yunus Can Gültekin
% Apr. 2021, Modified: Yunus Can Gültekin
% Feb. 2022, Modified: PAS with CCDM added. Kaiquan Wu
%--------------------------------------------------------------------------
% Input
% P: System Paramaters
% S: Signal Paramaters
%--------------------------------------------------------------------------
% Output
% P: System Paramaters - No updates
% S: Signal Paramaters
%    - .ShapAmps, the shaped amplitude sequence
%    - .FECinputBits, the shaped ampltude bits along with part of
%    informationbits for FEC encoding
%--------------------------------------------------------------------------

Sout = S;
Pout = P;

if strcmpi(P.SHA.Type,'ESS') || strcmpi(P.SHA.Type,'CCDM') || strcmpi(P.SHA.Type,'BESS')

    Nsym = P.Tx.Nsym;
    Rs   = P.SHA.Rs;
    Ks   = P.SHA.Ks;
    Nblk = P.SHA.Nblocks;

    dimNum = P.Sys.Npol*P.Sys.Nch;
    for rowind = 1:dimNum  % Loop over "streams" dimensions (polarizations and WDM channels)

        % For each row, the first 2*Nsym*Rs bits are to the amp.
        % shaper, the remaining 2*Nsym bits will be sign bits
        % Note that amplitudes are not 1, 3, 5..... but 0, 1, 2...
        if strcmpi(P.SHA.Type,'ccdm')     % CCDM
            shaped_ampIdx = zeros(Nblk*Pout.SHA.N,1);
            for blk = 1:Nblk
                pntBit = (blk-1)*Pout.SHA.Kshp+1:blk*Pout.SHA.Kshp;
                pntAmp = (blk-1)*Pout.SHA.N+1:blk*Pout.SHA.N;
                shaped_ampIdx(pntAmp) = CCDM_encode(Sout.InfoBits(rowind,pntBit).',Pout.SHA.N,Pout.SHA.Kshp,Pout.SHA.CC_Ni);
            end

        elseif strcmpi(P.SHA.Type,'ess') % ESS
            shaped_ampIdx = ESS_encode(P.SHA.Trellis,P.SHA.TrellisData,reshape(Sout.InfoBits(rowind,1:2*Nsym*Rs),Ks,P.SHA.Nblocks));
            shaped_ampIdx = reshape(shaped_ampIdx,[],1);

        elseif strcmpi(P.SHA.Type,'bess') % B-ESS
            shaped_ampIdx = ESS_BP_Shape(reshape(Sout.InfoBits(rowind,1:2*Nsym*Rs),Ks,P.SHA.Nblocks)',P.SHA.Trellis,P.SHA.N,P.SHA.L,P.SHA.na,P.SHA.Lmant,P.SHA.Kshp);
            shaped_ampIdx = (reshape(shaped_ampIdx,[],1)-1)/2;

        else
            error('Invalid PS type specification!');
        end

        % amplitude blocks
        TxIdx = shaped_ampIdx+sqrt(P.M)/2+1;
        Xsort = sort(P.P_1D.X);
        Sout.ShapAmps(rowind,:) = Xsort(TxIdx);

        % for PAS
        if ~strcmpi(P.FEC.CodeType,'none')
            % amplitude to bits
            ShapAmpBits = de2bi(Pout.P_1D.map(TxIdx),'left-msb');
            ShapAmpBits = reshape(ShapAmpBits(:,2:end).',[],Pout.FEC.Nblock);
            % the amplitude bits and info bits for FEC encoding
            pntBit0 = Pout.SHA.Nblocks*Pout.SHA.Kshp+1;
            ShapInfoBits = reshape(Sout.InfoBits(rowind,pntBit0:end).',[],Pout.FEC.Nblock);
            FECinputBits = [ShapAmpBits;ShapInfoBits];
            Sout.FECinputBits(rowind,1:Pout.FEC.k*Pout.FEC.Nblock) = reshape(FECinputBits,1,[]);
        end
    end


end

end

