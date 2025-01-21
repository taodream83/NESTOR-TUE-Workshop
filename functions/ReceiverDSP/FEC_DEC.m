function [Sout,Pout] = FEC_DEC(S,P)
%--------------------------------------------------------------------------
% This function implements FEC decoding.
% Jan. 2021, Created: Yunus Can Gültekin
% Feb. 2022, Modified: Kaiquan Wu. LDPC added
%--------------------------------------------------------------------------
% Input
% P: System Paramaters - No updates
%--------------------------------------------------------------------------
% S: Signal Paramaters
%    - .DecInfBitSeq, the decoded bit sequences
%--------------------------------------------------------------------------

Sout = S;
Pout = P;

switch lower(P.FEC.CodeType)

    % LDPC-coded transmission
    case {'dvbs2_ldpc','wifi_ldpc'}
        [Sout, Pout] = LDPCDecoder(S, P);
        Sout.RxSym.EstimatedBitsHD = Sout.DecInfBitSeq;

        % Product-coded transmission
    case 'pc'
        [Sout, Pout] = ProductDecoder(S, P);

end


end
