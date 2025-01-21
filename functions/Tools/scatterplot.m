function varargout = scatterplot(X,varargin)
% Improved scatterplot function compared to Matlab internal
% Supports both complex and real numbered inputs, in any matrix orientation.
% The optional second argument determines the marker size, defaults to 5 otherwise (Matlab default)
% If one of the dimensions of X is equal to 2 (complex numbers) or 4 (real numbers),
% dual graphs will be made, one for each polarization.
%
% Input:
% X             -Sequence or constellation to make a scatterplot of.
%                Both real and complex numbers are supported
% varagin       -Empty: defaults to markersize of 5
%               -Single number: Specifies markersize
%
% Output:
% varargout     -The scatterplot output figure handle
%
% Examples:
% scatterplot(Constellation); Plot constellation points with a default markersize of 5
% scatterplot(Constellation,10); Plot constellation points with a markersize of 10
%
% Sebastiaan Goossens
% Created: October 2023

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

markerSize = 5;
if nargin == 2
    if isnumeric(varargin{1}) && isscalar(varargin{1})
        markerSize = varargin{1};
    else
        warning('Invalid marker size input, defaulting to Matlab default (5)');
    end
end

h = figure;

if size(X,2) == 2 % Single polarization
    h.Position(3) = h.Position(4); % Make sure figure fits nicely around axes
    plot(X(:,1),X(:,2),'b.','MarkerSize',markerSize);

    % Adjust the limits (to mimic original scatterplot.m behavior)
    maxAll = max(max(X));
    limFact = 1.07;
    limits = maxAll*limFact;
    axis([-limits limits -limits limits]);

    title('Scatter plot');
    xlabel('In-Phase');
    ylabel('Quadrature');
    axis square;
    grid on;

elseif size(X,2) == 4 % Dual polarization
    h.Position(3) = 2*h.Position(4); % Make sure figure fits nicely around axes

    subplot(1,2,1);
    plot(X(:,1),X(:,2),'b.','MarkerSize',markerSize);

    % Adjust the limits (to mimic original scatterplot.m behavior)
    maxAll = max(max(X));
    limFact = 1.07;
    limits = maxAll*limFact;
    axis([-limits limits -limits limits]);
    sgtitle('Scatter plot');
    title('X-Polarization');
    xlabel('In-Phase');
    ylabel('Quadrature');
    axis square;
    grid on;

    subplot(1,2,2);
    plot(X(:,3),X(:,4),'b.','MarkerSize',markerSize);

    % Adjust the limits (to mimic original scatterplot.m behavior)
    axis([-limits limits -limits limits]);
    title('Y-Polarization');
    xlabel('In-Phase');
    ylabel('Quadrature');
    axis square;
    grid on;
end

if(nargout == 1)
    varargout(1) = {h};
end

end