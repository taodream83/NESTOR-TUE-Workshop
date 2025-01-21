function [Sout,Pout] = LDPCDecoder(S,P)
%% ------------------- LDPC Encoder --------------------%%
% A function that performs LDPC decoding and calculates the postBER FEC
% Kaiquan Wu, Feb. 2022
%--------------------------------------------------------------------------
% Input
% P: System Paramaters
% S: Signal Paramaters
%--------------------------------------------------------------------------
% Output
% P: System Paramaters
%    - .DecInfBitSeq, the decoded bit sequences
% S: Signal Paramaters - No updates
%--------------------------------------------------------------------------
Sout = S;
Pout = P;

pnt_dim = 1:P.Sys.Npol*length(P.Rx.CUT);
llr_dim = -S.RxSym.LLR(pnt_dim,:).';
for rowind = 1:size(llr_dim,2)  % Loop over "streams" dimensions (polarizations and WDM channels)
    llr = llr_dim(:,rowind);
    for nn=1:P.FEC.Nfecframe
        pnt_info    = (nn-1)*P.FEC.k+1:nn*P.FEC.k;
        pnt_coded   = (nn-1)*P.FEC.n+1:nn*P.FEC.n; 
        try
            Sout.DecInfBitSeq(rowind,pnt_info) = Pout.FEC.LDPCdec(llr(pnt_coded));
        catch
            Sout.DecInfBitSeq(rowind,pnt_info) = ldpcDecode(llr(pnt_coded),P.FEC.LDPCdeccfg,P.FEC.maxnumiter);
        end
    end
end
end

