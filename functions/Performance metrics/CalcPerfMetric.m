function R = CalcPerfMetric(R,S,P,varargin)

CUT=P.Rx.CUT;

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


for n = 1:length(varargin)
    %% Convert arguments to lowercase for case insensitive input
    switch lower(varargin{n})
        case 'snr'
            %% Calculate Signal-to-Noise Ratio (SNR)
            % Computes the effective SNR ssuming that the Tx
            % constellation is normalized to unit energy per polarization.

            % Author: Alex Alvarado
            % Modified by Gabriele Liga
            % Created: Jan 2020

            for chn = 1:length(TxChIdx)
                R.SNR(chn) = 1/mean(abs(S.SymSeq(TxChIdx(chn),:)-S.RxSym.EstimatedSymbols(chn,:)).^2); % SNR per 2D
            end

            R.SNR = reshape(R.SNR,P.Sys.Npol,length(P.Rx.CUT)); % Each column is a channel, and rows represent polarizations

            R.SNRdB = 10*log10(R.SNR);

        case 'mi'
            %% Calculate Mutual Information (MI)
            % The expression used is eq. (32) in "Achievable Information Rates for
            % Fiber Optics: Applications and Computations", by Alvarado et al., JLT
            % Jan. 2018.

            % Author: Sebastiaan Goossens
            % Created: February 2020

            if P.N == 2
                for chn = 1:length(TxChIdx)
                    X = S.SymSeq(TxChIdx(chn),:);
                    Y = S.RxSym.EstimatedSymbols(chn,:);

                    varZ = var(Y-X);
                    qYonX = (1/(pi*varZ))*exp(-abs(Y-X).^2/varZ);

                    qY = 0;
                    for ii=1:P.M
                        qY = qY + P.PrX(ii)*(1/(pi*varZ))*exp(-abs(Y-P.C(ii)).^2/varZ);
                    end

                    R.MI(chn) = (1/P.Tx.Nsym)*sum(log2(qYonX./qY));

                end
                if P.Sys.Npol == 2
                    R.MI = reshape(R.MI,2,length(P.Rx.CUT)); % Each column is a channel, and rows represent polarizations
                end

            elseif P.N == 4
                for chn = 1:length(TxChIdx)/2
                    X = S.SymSeq(TxChIdx(2*chn-1:2*chn),:);
                    Y = S.RxSym.EstimatedSymbols(2*chn-1:2*chn,:);

                    P_X = ones(1,P.M)/P.M;
                    varZ = mean(var(Y-X));
                    d = abs(Y(1,:)-X(1,:)).^2 + abs(Y(2,:)-X(2,:)).^2;
                    qYonX = (1/(pi*varZ))*exp(-d/varZ);

                    qY = 0;
                    for ii=1:P.M
                        d = abs(Y(1,:)-P.C(1,ii)).^2 + abs(Y(2,:)-P.C(2,ii)).^2;
                        qY = qY + P_X(ii)*(1/(pi*varZ))*exp(-d/varZ);
                    end

                    R.MI(chn) = (1/P.Tx.Nsym)*sum(log2(qYonX./qY));
                end
            end

        case 'gmi'
            %% Calculate Generalized Mutual Information (GMI)
            % Computes the GMI using LLRs calculated by the demapper.
            % The expression used is eq. (36) in "Achievable Information Rates for
            % Fiber Optics: Applications and Computations", by Alvarado et al., JLT
            % Jan. 2018.

            % Author: Bin Chen
            % Created: Jan 2020
            % Modified by Alex Alvarado
            % Modified by Gabriele Liga

            % Note that optimization over s is not actually included. It should be done
            % later. Small differences are probably observed.

            if P.N == 2
                for chn = 1:length(TxChIdx)
                    R.GMI(chn) = P.HX-1/P.Tx.Nsym*sum(log2(1+exp((-1).^(S.CodeIntBits(TxChIdx(chn),:)).*S.RxSym.LLR(chn,:))));
                end
                if P.Sys.Npol == 2
                    R.GMI = reshape(R.GMI,2,length(P.Rx.CUT)); % Each column is a channel, and rows represent polarizations
                end

            elseif P.N == 4
                for chn = 1:length(CUT)
                    R.GMI(chn) = P.m-1/P.Tx.Nsym*sum(log2(1+exp((-1).^(S.CodeIntBits(CUTidx(chn)+1,:)).*S.RxSym.LLR(chn,:))));
                end
            end

        case 'ser'
            %% Calculate Symbol Error Rate (SER)
            for chn = 1:length(TxChIdx)
                R.SER(chn) = sum(S.RxSym.EstimatedIndicesHD(chn,:) ~= S.TxIdx(TxChIdx(chn),:))/P.Tx.Nsym;
            end

            if P.Sys.Npol == 2
                R.SER = reshape(R.SER,2,length(P.Rx.CUT)); % Each column is a channel, and rows represent polarizations
            end

        case 'ber'
            %% Calculate Bit Error Rate (BER)
            if P.N == 2
                for chn = 1:length(TxChIdx)
                    if ~strcmpi(P.FEC.CodeType,'none')  % coded case
                        R.BER(chn) = sum(S.RxSym.EstimatedBitsHD(chn,:) ~= S.InfoBits(TxChIdx(chn),:))/(P.m*P.Tx.Nsym);
                    else
                        R.BER(chn) = sum(S.RxSym.EstimatedBitsHD(chn,:) ~= S.BitSeq(TxChIdx(chn),:))/(P.m*P.Tx.Nsym);
                    end
                end
                if P.Sys.Npol == 2
                    R.BER = reshape(R.BER,2,length(P.Rx.CUT)); % Each column is a channel, and rows represent polarizations
                end

            elseif P.N == 4
                for chn = 1:length(P.Rx.CUT)

                    if ~strcmpi(P.FEC.CodeType,'none') % coded case
                        R.BER(chn) = sum(S.RxSym.EstimatedBitsHD(chn,:) ~= S.InfoBits(chn,:))/(P.m*P.Tx.Nsym);
                    else
                        R.BER(chn) = sum(S.RxSym.EstimatedBitsHD(chn,:) ~= S.BitSeq(chn,:))/(P.m*P.Tx.Nsym);
                    end
                end
            end

        case 'q'
            %% Caclulate Q factor
            % Needs BER to be calculated first

            R.Q = 20*log10(sqrt(2)*erfcinv(2*R.BER)); % Q factor in dB

        otherwise
            %% Display a warning for unknown methods
            warning(['Unknown input: ' varargin{n}]);
    end
end

end
