function [Sout,Pout] = FEC_ENC(S,P)
%--------------------------------------------------------------------------
% This function implements FEC encoding.
% Jan. 2021, Created: Yunus Can Gültekin
% Apr. 2021, Modified: Yunus Can Gültekin
% Feb. 2022, Modified: Kaiquan Wu         -> LDPC Added
% Feb. 2022, Modified: Yunus Can Gültekin -> Simplified
%--------------------------------------------------------------------------
% Input
% P: System Paramaters
% S: Signal Paramaters
%--------------------------------------------------------------------------
% Output
% P: System Paramaters
%    - .Tx amd .FEC structs might be updated for product codes
% S: Signal Paramaters
%    - .CodedBits, the codeword sequences for each polarizations
%    - ...
%--------------------------------------------------------------------------

Sout = S;
Pout = P;

switch lower(P.FEC.CodeType)

    %========== LDPC-coded transmission
    case {'dvbs2_ldpc','wifi_ldpc'}
        [Sout, Pout] = LDPCEncoder(S,P);

        %========== Product-coded transmission
    case 'pc'
        [Sout, Pout] = ProductEncoder(S,P);

end

end

