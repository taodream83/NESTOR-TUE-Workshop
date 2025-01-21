function [Sout,Pout] = LDPCEncoder(S,P)
%--------------------------------------------------------------------------
% This function implements LDPC encoding.
% Feb. 2022. Kaiquan Wu
%--------------------------------------------------------------------------
% Input
% P: System Paramaters
% S: Signal Paramaters
%--------------------------------------------------------------------------
% Output
% P: System Paramaters. No updates
% S: Signal Paramaters
%    - .BitSeq, the coded bit sequence
%--------------------------------------------------------------------------

Sout = S;
Pout = P;

if ~strcmpi(P.SHA.Type,'none')
    inputBits = Sout.FECinputBits;
else
    inputBits = Sout.InfoBits;
end

dimNum = P.Sys.Npol*P.Sys.Nch;
for rowind = 1:dimNum  % Loop over "streams" dimensions (polarizations and WDM channels)
    for nn=1:P.FEC.Nfecframe
        pnt_info    = (nn-1)*P.FEC.k+1:nn*P.FEC.k;
        pnt_coded   = (nn-1)*P.FEC.n+1:nn*P.FEC.n;
        infbits     = inputBits(rowind,pnt_info).';
        block_coded = ldpcEncode(infbits, P.FEC.LDPCenccfg);
        Sout.BitSeq(rowind,pnt_coded) = block_coded;
    end
end

end