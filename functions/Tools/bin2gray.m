function [output,mapping] = bin2gray(x,modulation,order)
%BIN2GRAY Gray encode
%
%   Warning: BIN2GRAY will be removed in a future release. Use the
%   appropriate modulation object or function to remap constellation points
%   instead. See <a href="matlab:helpview('comm', 'CA_B2G')">this release note</a> for more information.
%
%   Y = BIN2GRAY(X,MODULATION,M) generates a Gray encoded output with the
%   same dimensions as its input parameter X.  The input X can be a scalar,
%   vector, matrix, or 3-D array.  MODULATION is the modulation type, which
%   can be a string equal to 'qam', 'pam', 'fsk', 'dpsk', or 'psk'.  M is
%   the modulation order that must be an integer power of two.
%
%   [Y,MAP] = BIN2GRAY(X,MODULATION,M) generates a Gray encoded output, Y and
%   returns its Gray encoded constellation map, MAP.  The constellation map is a
%   vector of numbers to be assigned to the constellation symbols.
%
%   If you are converting binary coded data to Gray coded data and modulating
%   the result immediately afterwards, you should use the appropriate modulation
%   functions with the 'gray' option, instead of BIN2GRAY.
%
%   EXAMPLE:
%     % To Gray encode a vector x with a 16-QAM Gray encoded constellation and
%     % return its map, use:
%     x=randi([0 15],1,100);
%     [y,map] = bin2gray(x,'qam',16);
%
%     % Obtain the symbols for 16-QAM
%     data = 0:15;
%     symbols = qammod(data,16,'bin');
%
%     % Plot the constellation
%     scatterplot(symbols);
%     set(get(gca,'Children'),'Marker','d','MarkerFaceColor','auto');
%     hold on;
%     % Label the constellation points according to the Gray mapping
%     for jj=1:16
%       text(real(symbols(jj))-0.15,imag(symbols(jj))+0.15,...
%       dec2base(map(jj),2,4));
%     end
%     set(gca,'yTick',(-4:2:4),'xTick',(-4:2:4),...
%      'XLim',[-4 4],'YLim',...
%      [-4 4],'Box','on','YGrid','on', 'XGrid','on');
%
%   See also GRAY2BIN, PSKMOD, QAMMOD, PAMMOD, FSKMOD.

%   Copyright 1996-2021 The MathWorks, Inc.

%#codegen

%% Begin validating inputs

% Typical error checking.
narginchk(3, 3)

%Validate numeric x data
if isempty(x)
    coder.internal.error('comm:bin2gray:InputEmpty');
end

% x must be a scalar, vector, matrix or 3D array
if length(size(x)) > 3
    coder.internal.error('comm:bin2gray:InputDimensions');
end

% x must be a finite non-negative integer
xVec = x(:);
if (any(xVec<0) || any(isinf(xVec)) || (~isreal(x)) || any(floor(xVec) ~= xVec))
    coder.internal.error('comm:bin2gray:InputError');
end

% Validate modulation type
if (~comm.internal.utilities.isCharOrStringScalar(modulation)) || (~strcmpi(modulation,'QAM')) && (~strcmpi(modulation,'PSK'))...
        && (~strcmpi(modulation,'FSK')) && (~strcmpi(modulation,'PAM')) && (~strcmpi(modulation,'DPSK'))
    coder.internal.error('comm:bin2gray:ModulationTypeError');
end

%Validate modulation order
if (order < 2) || (isinf(order) || ...
        (~isreal(order)) || (floor(log2(order)) ~= log2(order)))
    coder.internal.error('comm:bin2gray:ModulationOrderError');
end

% Check for overflows - when x is greater than the modulation order
if (max(xVec) >= order)
    coder.internal.error('comm:bin2gray:XError');
end

%% Start Gray code conversion
[output, mapping] = bin2grayInternal(x, modulation, order);
end

function [output, mapping] = bin2grayInternal(x, modulation, order)
% COMM.INTERNAL.UTILITIES.BIN2GRAY The algorithm used by BIN2GRAY function
%
%   Y = COMM.INTERNAL.UTILITIES.BIN2GRAY(X, MODULATION, M)
%   [Y,MAP] = COMM.INTERNAL.UTILITIES.BIN2GRAY(X, MODULATION, M)
%
%   The input and output arguments are as listed in BIN2GRAY function.
%
%   Assumptions & Notes:
%   1) Supports codegen
%   2) Supports fixed-point input-output only for QAM modulation. The
%   fixed-point datatype of Y and MAP is same as that of input X.
%   3) Supports scalar, vector, matrix and 3 dimensional array input.
%   4) This is an internal function which does not do any input
%   validations. It is users responsibility to ensure that the above
%   conditions are met and that the input arguments are valid as per the
%   specifications of BIN2GRAY function.

%   Copyright 2015-2018 The MathWorks, Inc.

%#codegen

switch lower(modulation)

    case {'psk','pam','fsk','dpsk'}

        % Calculate Gray table
        j = int32((0:order-1)');
        mapping = cast(bitxor(j,bitshift(j,-1)),'like',x);

        % Format output and translate x (map) i.e. convert to Gray
        if isrow(x)
            % Assure that the output, if one dimensional,
            % has the correct orientation
            tmpMapping = mapping';
            output = tmpMapping(x+1);
        else
            output = mapping(x+1);
        end

    case {'qam'}

        [mapping, symbolIndex] = getGrayMapping(order, x);

        % Format output and translate x (map) i.e. convert to Gray
        if isrow(x)
            % symbolIndex AND x, both, are vectors AND they have
            % different orientations. As symbolIndex is always a column
            % vector, this degenerates into isrow(x) check.
            tmpSymbolIndex = symbolIndex';
            output = tmpSymbolIndex(x + cast(1,'like',x));
        else
            output = symbolIndex(x + cast(1,'like',x));
        end

    otherwise
        coder.internal.error('comm:bin2gray:ModulationTypeUnknown');
end
end

function [mapping, outSymbolIndex] = getGrayMapping(M, protoType)
% COMM.INTERNAL.QAM.GETGRAYMAPPING Generate Gray code for QAM
%
%   MAP = COMM.INTERNAL.QAM.GETGRAYMAPPING(M, PROTOTYPE)
%   M is modulation order and must be an integer power of 2.
%   The output has same datatype (and complexity) as PROTOTYPE.
%   MAP is the Gray encoded constellation map, as described by bin2gray
%   and gray2bin functions.
%
%   [MAP, SYMBOLORDER] = COMM.INTERNAL.QAM.GETGRAYMAPPING(M, PROTOTYPE)
%   SYMBOLORDER is the intermediate symbol order that is used by
%   bin2gray.
%
%   Assumptions & Notes:
%   1) Supports codegen.
%   2) Supports fixed-point input-output.
%   3) The datatype of MAP & SYMBOLORDER is same as that of input
%   PROTOTYPE. For efficiency, pass in PROTOTYPE as a scalar value of
%   the right type.
%   4) MAP & SYMBOLORDER are column vectors with M unique integer
%   valued elements in the range [0, M-1].
%   4) This is an internal function which does not do any input
%   validations. It is users responsibility to ensure that the above
%   conditions are met and that the input arguments are valid.

%   Copyright 2015-2021 The MathWorks, Inc.

%#codegen

% Number of bits per symbol
nBits = int32(0);
if isfi(protoType)
    tmp = int32(M);
    while tmp ~= 1
        nBits(1) = nBits(1) + 1;
        tmp = bitshift(tmp, -1);
    end
    binMapping = cast((0:M-1)', 'like', protoType);
    SymbolIndex = zeros(size(binMapping), 'like', binMapping);
else
    nBits(1) = log2(M);
    binMapping = int32((0:M-1)');
    SymbolIndex = zeros(size(binMapping), 'int32');
end

% mapping = zeros(M, 1, 'like', protoType);
mapping = cast(zeros(M, 1), 'like', protoType);
outSymbolIndex = zeros(M, 1, 'like', mapping);
tmpM = cast(M, 'like', binMapping);
tmpOne = cast(1,'like',tmpM);

if bitget(nBits,1)
    % Cross constellation

    nBitsI = bitshift(nBits+1, -1); % divide by 2
    nBitsQ = bitshift(nBits-1, -1);

    symbolI = bitshift(binMapping, -nBitsQ);
    symbolQ = bitand(binMapping, cast(bitshift(tmpM-tmpOne,-nBitsI),'like',binMapping));

    i = 1;
    while i < nBitsI
        tmpI = symbolI;
        tmpI = bitshift(tmpI,-i);
        symbolI = coder.sameSizeBinaryOp(@bitxor,symbolI,tmpI);
        % i takes on values 1,2,4,8,...,2^n - n is an integer
        i = i + i;
    end

    i = 1;
    while i < nBitsQ
        tmpQ = symbolQ;
        tmpQ = bitshift(tmpQ,-i);
        symbolQ = coder.sameSizeBinaryOp(@bitxor, symbolQ, tmpQ);
        % i takes on values 1,2,4,8,...,2^n - n is an integer
        i = i + i;
    end

    SymbolIndex(:) = bitshift(symbolI,nBitsQ) + symbolQ;

else % square constellation

    nBitsBy2 = bitshift(nBits,-1);
    symbolI = bitshift(binMapping, -nBitsBy2);
    symbolQ = bitand(binMapping, cast(bitshift(tmpM-tmpOne,-nBitsBy2),'like',binMapping));

    i = 1;
    while i < nBitsBy2
        tmpI = symbolI;
        tmpI = bitshift(tmpI,-i);
        symbolI = coder.sameSizeBinaryOp(@bitxor,symbolI,tmpI);

        tmpQ = symbolQ;
        tmpQ = bitshift(tmpQ,-i);
        symbolQ = coder.sameSizeBinaryOp(@bitxor, symbolQ, tmpQ);

        % i takes on values 1,2,4,8,...,2^n - n is an integer
        i = i + i;
    end

    SymbolIndex(:) = bitshift(symbolI, nBitsBy2) + symbolQ;

end

outSymbolIndex(:) = SymbolIndex;

% Make sure that mapping is a vector, when used to name the symbols
% column-wise starting from left upper corner, results in a gray mapped
% constellation.
mapping(outSymbolIndex + cast(1,'like',outSymbolIndex)) = 0:M-1;

end