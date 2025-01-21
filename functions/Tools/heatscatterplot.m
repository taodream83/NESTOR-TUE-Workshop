function varargout = heatscatterplot(X,varargin)
% Scatterplot with heatmap functionality. Uses 2d histogram internally to determine weight of the heatmap.
% Supports both complex and real numbered inputs, in any matrix orientation.
% If no second argument is given, the algorithm will automatically recognize if the input is on a grid up
% to and including 4096QAM, and bin it accordingly.
% Manual can be done by specifying a number as the second argument, which
% will be the number of bins. Another method is to specify the constellation.
% If one of the dimensions of X is equal to 2 (complex numbers) or 4 (real numbers),
% dual graphs will be made, one for each polarization.
%
% Input:
% X             -Sequence or constellation to make a heatscatterplot of.
%                Both real and complex numbers are supported
% varagin       -Empty: input provides scatterplot with predefined number of bins
%               -Single number: specifies number of bins
%               -Array: specifies probabilities (only use with constellation plot)
%
% Output:
% varargout     -The heatscatterplot output figure handle
%
% Examples:
% heatscatterplot(TxSymbols);      Plot transmitted symbols with grid recognition
% heatscatterplot(TxSymbols,16);   Plot transmitted symbols with grid recognition and a grid of 16 bins per axis
% heatscatterplot(RxSymbols);      Plot received symbols with automatic predetermined binning
% heatscatterplot(RxSymbols,1000); Plot received symbols with automatic binning on 1000 bins per axis
% heatscatterplot(Constellation,Probabilities); Plot constellation with bin location determined by the constellation and bin hight by the probabilities
%
% Sebastiaan Goossens
% Created: November 2018
% Updated: October 2023

add = isreal(X)*2; % Double amount of dimensions if real numbers are provided
if size(X,1) > 2+add && size(X,2) > 2+add % Treat as a single polarization matrix
    X = X(:);
end
if size(X,2) > 2+add % Force column vectors otherwise
    X = X.';
end

if ~isreal(X)
    if size(X,2) == 2
        X = [real(X(:,1)) imag(X(:,1)) real(X(:,2)) imag(X(:,2))];
    else
        X = [real(X) imag(X)];
    end
end

% Use up to idxNum elements to determine grid
idxNum = 10000;
if numel(X) < idxNum
    lastIdx = numel(X);
else
    lastIdx = idxNum;
end

% Select amplitudes of up to idxNum elements to not waste too much time on this
Z = reshape(X(1:lastIdx),[],1);

% Calculate minimum distance between constellation points
minDiff = min(nonzeros(diff(sort(Z))));
Z = abs(Z);
% Calculate number of unique amplitudes
Znum = numel(unique(Z));
Zmax = max(abs(X(:))); % Use X for determining max, not Z

% Set tolerance for binning with floating point numbers (single precision epsilon is 1.1921e-07)
tol = 1e-6;

% Calculate if all points are on the grid (with some tolerance for floating point imprecision)
onGrid = ~any(mod(Z,minDiff/2) > tol) && minDiff > tol;

% If only the symbol sequence is given as input:
if nargin == 1
    % If on the grid and limited amount of amplitudes (allows up to 4096 QAM)
    % Else a noisy signal is received, set default number of bins
    if onGrid && Znum <= 32
        edges = -(Zmax+minDiff/2):minDiff:((Zmax+minDiff/2)+tol); % +tol for floating point imprecision
    elseif Znum <= 32
        nBins = 44; % Largest value for which histogram2 paints black borders around bins
        edges = linspace(-Zmax-tol,Zmax+tol,nBins+1);
    else
        nBins = 250; % Default value, looks nice enough in most cases
        edges = linspace(-Zmax-tol,Zmax+tol,nBins+1);
    end
else % either a bin number or probabilities are given
    if numel(varargin{1}) == 1 % Assume this means a number of bins is provided
        nBins = varargin{1};
        if onGrid && nBins <= 64 % Assume that a nicely spaced constellation grid is provided (e.g. QAM)
            edges = -(nBins/2*minDiff):minDiff:(nBins/2*minDiff)+tol; % +tol for floating point imprecision
        else % Use provided number of bins
            edges = linspace(-Zmax-tol,Zmax+tol,nBins+1);
        end
    else % Assume the input is probabilities corresponding to the constellation
        if ~(onGrid && size(X,2) == 2)
            error('Predefined probabilities currently only possible in a square, single-pol constellation');
        else
            P = varargin{1};
            if length(X) ~= numel(P)
                error('If the second argument is a probability array, the first argument must be the constellation (of equal length).');
            end
            edges = -(sqrt(length(X))/2*minDiff):minDiff:(sqrt(length(X))/2*minDiff)+tol; % +tol for floating point imprecision
        end
    end
end

h = figure;

if size(X,2) == 2 % Single polarization
    h.Position(3) = h.Position(4); % Make sure figure fits nicely around axes

    if exist('P','var')
        [~,idx] = sort(complex(X(:,2),X(:,1)),'ComparisonMethod','real'); % Cannot sort by imag first unfortunately, so use complex() workaround
        P = P(idx); % Sort probabilities for use with histogram2()
        histogram2('XBinEdges',edges,'YBinEdges',edges,'BinCounts',reshape(P,sqrt(length(X)),sqrt(length(X))),'FaceColor','flat');
    else
        histogram2(X(:,1),X(:,2),edges,edges,'Normalization','probability','FaceColor','flat');
        if ~onGrid
            % Adjust the limits (to mimic original scatterplot.m behavior)
            maxAll = max(max(X));
            limFact = 1.07;
            limits = maxAll*limFact;
            axis([-limits limits -limits limits]);
        end
    end
    title('Scatter plot');
    xlabel('In-Phase');
    ylabel('Quadrature');
    axis square;
    grid on;
    view(2);

elseif size(X,2) == 4 % Dual polarization
    h.Position(3) = 2*h.Position(4); % Make sure figure fits nicely around axes

    subplot(1,2,1);
    histogram2(X(:,1),X(:,2),edges,edges,'FaceColor','flat');
    if ~onGrid
        % Adjust the limits (to mimic original scatterplot.m behavior)
        maxAll = max(max(X));
        limFact = 1.07;
        limits = maxAll*limFact;
        axis([-limits limits -limits limits]);
    end
    sgtitle('Scatter plot');
    title('X-Polarization');
    xlabel('In-Phase');
    ylabel('Quadrature');
    axis square;
    grid on;
    view(2);

    subplot(1,2,2);
    histogram2(X(:,3),X(:,4),edges,edges,'FaceColor','flat');

    if ~onGrid
        % Adjust the limits (to mimic original scatterplot.m behavior)
        axis([-limits limits -limits limits]);
    end
    title('Y-Polarization');
    xlabel('In-Phase');
    ylabel('Quadrature');
    axis square;
    grid on;
    view(2);
end

if(nargout == 1)
    varargout(1) = {h};
end

end