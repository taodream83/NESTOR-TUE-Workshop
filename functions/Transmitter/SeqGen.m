function [varargout] = SeqGen(S,P)
%--------------------------------------------------------------------------
% This function generates trasmitted symbol sequences.
% Jan. 2019, Bin Chen
% Jan. 2021, Yunus Can Gültekin
% Apr. 2021, Yunus Can Gültekin
% Feb. 2022, Kaiquan Wu. PAS enabled
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

Sout = S;
Pout = P;
Nsym = P.Tx.Nsym;

if P.N == 2
    %% Unshaped & uncoded transmission
    if strcmpi(P.FEC.CodeType,'none') && strcmpi(P.SHA.Type,'none')
        CodeIntBits = Sout.BitSeq; % no interleaving
        % Indices for mapper
        pad = mod(size(CodeIntBits,2),P.m);
        if ~pad == 0
            CodeIntBits = [CodeIntBits,zeros(P.Sys.Npol*P.Sys.Nch,pad)]; % Zero-padding to map to the closest number of symbols containing the coded bit sequence
        end
        TxIdx = reshape(b2d(reshape(CodeIntBits.',P.m,[]).'),Nsym,[]).'+1;

        %% Unshaped & coded transmission
    elseif ~strcmpi(P.FEC.CodeType,'none') && strcmpi(P.SHA.Type,'none')
        switch lower(P.FEC.CodeType)
            case {'dvbs2_ldpc','wifi_ldpc'} %LDPC
                CodeIntBits = Sout.BitSeq; % no interleaving
            case 'pc' %BCH/SCC or other HD-FEC
                Nsym = Pout.Tx.Nsym;
                CodeIntBits = Sout.CodedBits; % no interleaving
        end
        % Padding
        Pout.FEC.pad = mod(size(CodeIntBits,2),P.m);
        if Pout.FEC.pad ~= 0
            CodeIntBits = [CodeIntBits,zeros(P.Sys.Npol*P.Sys.Nch,Pout.FEC.pad)]; % Zero-padding to map to the closest number of symbols containing the coded bit sequence
        end
        % Indices for mapper
        TxIdx = reshape(b2d(reshape(CodeIntBits.',P.m,[]).'),Nsym,[]).'+1;

        %% Shaped & uncoded transmission
    elseif strcmpi(P.FEC.CodeType,'none') && ~strcmpi(P.SHA.Type,'none')
        Nblk = Pout.SHA.Nblocks;
        Rs   = P.SHA.Rs;

        dimNum = P.Sys.Npol*P.Sys.Nch;
        for rowind = 1:dimNum % Loop over "streams" dimensions (polarizations and WDM channels)
            shap_amps = Sout.ShapAmps(rowind,1:Nblk*P.SHA.N);
            shap_ampIdx = (round(shap_amps./min(abs(shap_amps)))+1)/2-1; % convert amplitudes to index
            shap_amps_g = bin2gray(shap_ampIdx,'pam',P.SHA.na);
            % With the following 2 lines, the first half of the
            % amps/signs are put in I and the second half are in Q.
            % In the coded case, this would lead to coded bits being
            % separated many chan. uses, i.e., amp.-level interleaving.
            syms_quad = sqrt(P.M/4)*shap_amps_g(1:numel(shap_amps_g)/2) + shap_amps_g(numel(shap_amps_g)/2+1:end); % Symbols for single quadrant

            TxIdx(rowind,:) = 2*(P.M/4)*Sout.InfoBits(rowind,2*Nsym*Rs+1:2*Nsym*Rs+Nsym)...
                + (P.M/4)*Sout.InfoBits(rowind,2*Nsym*Rs+Nsym+1:end) + syms_quad + 1;
            CodeIntBits(rowind,:) = reshape(d2b(TxIdx(rowind,:)-1, P.m)',[],1)';
        end

        %% Shaped & coded transmission, probabilistic amplitude shaping (PAS)
    elseif ~strcmpi(P.FEC.CodeType,'none') && ~strcmpi(P.SHA.Type,'none')
        dimNum = P.Sys.Npol*P.Sys.Nch;
        for rowind = 1:dimNum % Loop over "streams" dimensions (polarizations and WDM channels)
            % extract sign bits
            CodeIntBits = reshape(Sout.BitSeq(rowind,:),[],Pout.FEC.Nblock);
            parityBits = CodeIntBits(1+Pout.FEC.k:end,:);
            pntBit0 = Pout.SHA.Nblocks*Pout.SHA.Kshp+1;
            ShapInfoBits = reshape(Sout.InfoBits(rowind,pntBit0:end).',[],Pout.FEC.Nblock);
            signBits = [parityBits;ShapInfoBits];
            % combine sign bits and shaped amplitudes to be shaped 2D rectangular QAM
            amps = reshape(Sout.ShapAmps(rowind,:),[],Pout.FEC.Nblock);
            symsPAM = reshape(amps.*(-1).^(1-signBits),2,[]);
            Sout.SymSeq(rowind,:) = symsPAM(1,:)+symsPAM(2,:).*1i;
        end
        Sout.CodeIntBits = Sout.BitSeq;
        varargout{1} = Sout;
        return
    end

elseif P.N == 4
    if strcmpi(P.FEC.CodeType,'none') && strcmpi(P.SHA.Type,'none')
        CodeIntBits = Sout.BitSeq; % no interleaving
        % Indices for mapper
        pad = mod(size(CodeIntBits,2),P.m);
        if ~pad == 0
            CodeIntBits = [CodeIntBits,zeros(P.Sys.Npol*P.Sys.Nch,pad)]; % Zero-padding to map to the closest number of symbols containing the coded bit sequence
        end
        TxIdx = reshape(b2d(reshape(CodeIntBits.',P.m,[]).'),Nsym,[]).'+1;
    elseif strcmpi(P.FEC.CodeType,'none') && strcmpi(P.SHA.Type,'ideal')
        TxIdx = sum(rand(1,Nsym) > Pout.P_apriori_cumulative) + 1;
        CodeIntBits = reshape(d2b(TxIdx-1, P.m).',1,[]);
    end
end

%% Integer-to-symbol mapping
Sout.TxIdx = TxIdx;
if P.N == 2
    Sout.SymSeq = P.C(Sout.TxIdx);
elseif P.N == 4
    Sout.SymSeq = zeros(P.Sys.Nch*P.Sys.Npol,Nsym);
    for n = 1:P.Sys.Nch
        Sout.SymSeq(n*2-1:n*2,:) = P.C(:,TxIdx(n,:));
    end
end

if ~isfield(P,'FEC') || strcmpi(P.FEC.CodeType,'none') % no FEC
    Sout.CodeIntBits = CodeIntBits;
    varargout{1} = Sout;
else
    Sout.CodeIntBits = CodeIntBits;
    varargout{1} = Sout;
    varargout{2} = Pout;
end

end


