function Pout = ModGen(P)
% This function returns the constellation points for M-PSK and M-PAM for simple testing
% The format for the constellation is an N by M matrix (double)
% The format for the labeling is an M by m matrix (double)
% Author: Bin Chen
% Created: Jan. 2019
% Modified:
% Gabriele Liga, Oct. 2019
% Alex Alvarado, Dec. 2019
% Yunus an Gültekin, Apr. 2021
% Kaiquan Wu, Feb. 2022. Enable CCDM in the case of PAS added
%--------------------------------------------------------------------------
% Input
% P: System Paramaters
%--------------------------------------------------------------------------
% Output
% P: System Paramaters
%    - .X, the constellation set corresponding to each dimension (inphase/quadrature, X/Y pol)
%    - .C, 2D/4D Complex-valued constellation
%    - .map, the mappping rule from constellation X to binary labeling (the labeling are shown in decimal)
%    - .M, modulation cardinality and order
%    - .m, modulation order
%    - .N, the constellation dimension
%    - .PrX, probability of constellation,
%    - .HX, entropy of constellation,
%    - .P_1D, the signal struct for 1D PAM modulation when we do PAS for 2D
%--------------------------------------------------------------------------

% in case the shaping field is not defined
if ~isfield(P,'SHA')
    P.SHA.Type = 'none';
end

Pout = P;
cons = P.ModFormat;
M = P.M ;
m = log2(M);
Pout.m = m;
Pout.PrX = ones(M,1)/M; % Gets overwritten with PS
Pout.HX = m; % Gets overwritten with PS

%% 2D Modulation
if P.N == 2
    switch cons %case {'PAM','PSK','QAM'}. Labeling is the 'BRGC'
        case'PAM'
            Delta=sqrt(3/(M^2-1));
            X=[Delta*(-(M-1):2:(M-1));zeros(1,M)];
            if isfield(P,'BinLabeling') && ~isequal(P.BinLabeling,'Gray')
                %% To be completed with well-known labelings alternative to Gray
            else
                [~,map]=bin2gray(0:M-1,'pam',M);
            end
        case'PSK'
            X=[cos((1:M)*2*pi/M-pi/M);sin((1:M)*2*pi/M-pi/M)];
            if isfield(P,'BinLabeling') && ~isequal(P.BinLabeling,'Gray')
                %% To be completed with well-known labelings alternative to Gray
            else
                [~,map]=bin2gray(0:M-1,'psk',M);
            end

        case 'QAM'
            if strcmpi(P.SHA.Type,'none') || strcmpi(P.SHA.Type,'ideal')
                Delta=sqrt(3/(M-1));
                x=Delta*(-(sqrt(M)-1):2:(sqrt(M)-1));
                I=repmat(x,sqrt(M),1);
                Q=I.';
                X=zeros(2,M);
                for nn=1:M
                    X(1,nn)=I(nn);
                    X(2,nn)=-Q(nn);    % Assign position of each constellation point in the array following graphical order (linear index spanning across columns)
                end
                if isfield(P,'BinLabeling') && ~isequal(P.BinLabeling,'Gray')
                    % To be completed with well-known labelings alternative to Gray
                else
                    [~,map]=bin2gray(0:M-1,'qam',M);
                end
            else
                map = 0:P.M-1;
                % PAS uses PAM for each dimension, get the 1D PAM
                % modulation sturct P_1D, which will be used later for
                % sequence generation and demapper
                P_1D = P;
                P_1D.ModFormat = 'PAM';P_1D.M=sqrt(P_1D.M); P_1D.m=log2(P_1D.M);P_1D.SHA.Type = 'none';
                P_1D = ModGen(P_1D);
                P_1D.X = P_1D.X(1,:); % WORKAROUND: Force P_1D.X to be 1D

                if strcmpi(P.SHA.Type,'ess')
                    if ~isfield(P.SHA,'L')
                        Pout.PS.L = ESS_calc_rate(P.SHA.Rs,P.SHA.N,P.M_dim,P.SHA.Lmant);
                    end
                    [Pout.SHA.Trellis, Pout.SHA.TrellisData, Pout.SHA.Kshp, Pout.SHA.PrA_1D] = ESS_generateTrellis(P.SHA.N,P.SHA.L,P.SHA.na,P.SHA.Lmant);

                elseif strcmpi(P.SHA.Type,'bess')
                    [Pout.SHA.Trellis, Pout.SHA.Kshp] = BESS_generateTrellis (P.SHA.N,P.SHA.L,P.SHA.na,P.SHA.Lmant,P.SHA.hi,P.SHA.s,P.SHA.wi);
                    if Pout.SHA.Kshp < P.SHA.Ks
                        error('Parameters of the trellies do not satisfy the target rate.');
                    else
                    end
                    [Pout.SHA.PrA_1D,~,~] = BESS_estimatePMF(P.SHA.N,P.SHA.L,P.SHA.na,P.SHA.Lmant,P.SHA.hi,P.SHA.s,P.SHA.wi);
                    P_1D.SHA.PrA_1D = Pout.SHA.PrA_1D';
                    
                end
               
                % the variables X, and map for 2D will not be used for
                % sequence generation and demapper
                [X, Pout.PrX, Pout.HX] = PS_generateConst(Pout);
            end
        case 'custom'
            X = P.CustomConstellation;
            map = P.CustomMapping;
        otherwise
            error('Invalid constellation type');
    end

    %% 4D Modulation
elseif P.N == 4
    if strcmpi(cons,'custom')
        X = P.CustomConstellation;
        map = P.CustomMapping;
    else
        M_2D = sqrt(M); % Define amount of points in 2D
        switch cons
            case 'PAM' % Create 4D equivalent of 2D-PAM
                C_2D = pammod(0:(M_2D-1),M_2D,0,'gray'); % Define C as standard M_2D-PAM
            case 'PSK' % Create 4D equivalent of 2D-PSK
                C_2D = pskmod(0:(M_2D-1),M_2D,pi/M,'gray'); % Define C as standard M_2D-PSK
            case 'QAM' % Create 4D equivalent of 2D-QAM
                C_2D = qammod(0:(M_2D-1),M_2D,'gray'); % Define C as standard M_2D-QAM
            otherwise
                error('Invalid constellation type');
        end
        C_4D = [repelem(C_2D,M_2D);repmat(C_2D,1,M_2D)]; % Define as 4D constellation
        X = [real(C_4D(1,:));imag(C_4D(1,:));real(C_4D(2,:));imag(C_4D(2,:));]; % Convert to separate dimensions
        map = 1:P.M; % Constellation is already Gray labeled, set linear order
    end
end

Pout.N=size(X,1);
Pout.map=map;

%% Probabilistic Shaping logic
% TODO: Rework all PS logic so it is completely separated from the
% constellation geometry definitions
if strcmpi(P.SHA.Type,'ideal')
    if strcmpi(P.SHA.distribution,'custom')
        Pout.PrX = P.SHA.CustomProbabilities;
    end
    Pout.PrX = Pout.PrX(:)/sum(Pout.PrX(:)); % Make sure this is normalized and in column vector
    Pout.PrX_cumulative = cumsum(Pout.PrX);
end

%% Permute constellation points to match decimal labeling order
[~,idx_sort] = sort(map);
try
    X = X(:,idx_sort);
catch
    X = X(idx_sort);
end
Pout.X = normalize_Es(X,Pout.PrX,Pout.N); % unit energy constellation

if ~strcmpi(P.SHA.Type,'none') && ~strcmpi(P.SHA.Type,'ideal')
    % If we have PS, align the 1D PAM set as the normalized shaped 2D QAM
    % projected on one dimension
    PAM_C = sort(unique(Pout.X(1,:)));
    [~,idx_sort] = sort(P_1D.X);
    [~,r] = sort(idx_sort);
    P_1D.X = PAM_C(r);
    Pout.P_1D = P_1D;
end

%% 2D/4D real to 1D/2D complex mapping
if Pout.N == 2
    Pout.C = (Pout.X(1,:) + 1i*Pout.X(2,:));               % 2D Complex-valued constellation
elseif Pout.N == 4
    Pout.C = (Pout.X([1 3],:) + 1i*Pout.X([2 4],:));       % 4D Complex-valued constellation
end

%% Find Subconstellations defined by the labeling (use for demapper)
Pout.Ik0 = zeros(P.M/2,Pout.m);
Pout.Ik1 = zeros(P.M/2,Pout.m);

Lbin=de2bi(0:P.M-1,'left-msb');

for kk=1:Pout.m
    pntr0=1;
    pntr1=1;
    for i=1:P.M
        if Lbin(i,kk)==0
            Pout.Ik0(pntr0,kk)=i;
            pntr0=pntr0+1;
        else
            Pout.Ik1(pntr1,kk)=i;
            pntr1=pntr1+1;
        end
    end
end
end

%% Functions
function X = normalize_Es(X,PX,N)
Es=sum(X.^2,1)*PX(:);
X=X./sqrt(Es);
if N == 4 % To have unity enery per polarization for 4D constellations
    X=X.*sqrt(2);
end
end

function [X, PB, HX] = PS_generateConst(P)
% Generate QAM constellation and probabilities from 1D parameters
if P.SHA.na ~= 1
    [x,idx] = bin2gray(0:P.SHA.na-1,'pam',P.SHA.na);
    [~,idx] = sort(idx);
    idx = idx - 1;
else
    x = 0;
    idx = 0;
end

x = x + P.SHA.na.*rot90(x);            % Convert to 2D plane to get gray coded quadrant
x = [fliplr(P.M/4+x) x];       % Flip vertically and concatenate to generate upper half of symbolmap
x = [x;flipud(2*P.M/4+x)];     % Flip horizontally and concatenate to generate full symbolmap
Xc = qammod((0:P.M-1)',P.M,x(:));
X(1,:) = real(Xc);
X(2,:) = imag(Xc);
PB = (P.SHA.PrA_1D(idx+1)).*(P.SHA.PrA_1D(idx+1).')/4;
PB = repmat(PB(:),4,1);
HX = -sum(PB(PB > 0).*log2(PB(PB > 0)));
end
