function Sout = ConstSync(S,P)
%CONSTSYNC
% Sout = ConstSync(S,P) Recovers the original phase and scaling of a signal
% The following methods are available: 'xcorr'
%
% Input:
% S     -Signal structure
% P     -Parameter structure
%
% Output:
% Sout  - Output signal structure

% Author: different authors per method

switch P.Rx.ConstSyncMethod

    case 'xcorr'
        %% XCORR
        % Returns the normalized circular covariance of 2 complex sequences and
        % calculates the delay between 2 sequences
        %
        % INPUTS:
        % Input sequences can be 2xN matrices in which case the
        % function calculates the correlation for each row of the matrix
        %
        % OUTPUTS:
        % C   - complex correlation
        % lag - delays
        %
        % Author: Gabriele Liga, April 2019
        % Author: Astrid Barreiro, January 2020


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


        for chn = 1:length(TxChIdx)

            x = S.SymSeq(TxChIdx(chn),:);
            y = S.RxSym.EstimatedSymbols(chn,:);
            C = ifft(fft(x,[],2).*conj(fft(y,[],2)),[],2);

            lag = zeros(size(C,1),1);
            for ii = 1:size(C,1)
                lag(ii) = find(abs(C(ii,:)) == max(abs(C(ii,:)))) - 1;
            end

            S.RxSym.EstimatedSymbols(chn,:) = circshift(S.RxSym.EstimatedSymbols(chn,:),[0,lag])*1/(abs(C(lag+1))/length(C))*exp(1i*angle(C(lag+1)));

        end

    otherwise % If no correct parameter is given
        error("ConstSync: Please input one of the following as parameter: 'xcorr'")
end

%% Set output variables
Sout = S;
end