function y = dbmean(x,varargin)
%DBMEAN   Average or mean value of an array containing logarighmic powers.
%   This documentation is taken from the Matlab implementation of MEAN.
%   S = DBMEAN(X) is the mean value of the elements in X if X is a vector.
%   For matrices, S is a row vector containing the mean value of each
%   column.
%   For N-D arrays, S is the mean value of the elements along the first
%   array dimension whose size does not equal 1.
%
%   DBMEAN(X,'all') is the mean of all elements in X.
%
%   DBMEAN(X,DIM) takes the mean along the dimension DIM of X.
%
%   DBMEAN(X,VECDIM) operates on the dimensions specified in the vector
%   VECDIM. For example, DBMEAN(X,[1 2]) operates on the elements contained
%   in the first and second dimensions of X.
%
%   S = DBMEAN(...,TYPE) specifies the type in which the mean is performed,
%   and the type of S. Available options are:
%
%   'double'    -  S has class double for any input X
%   'native'    -  S has the same class as X
%   'default'   -  If X is floating point, that is double or single,
%                  S has the same class as X. If X is not floating point,
%                  S has class double.
%
%   S = DBMEAN(...,NANFLAG) specifies how NaN (Not-A-Number) values are
%   treated. The default is 'includenan':
%
%   'includenan' - the mean of a vector containing NaN values is also NaN.
%   'omitnan'    - the mean of a vector containing NaN values is the mean
%                  of all its non-NaN elements. If all elements are NaN,
%                  the result is NaN.
%
%   Example:
%       X = [1 2 3; 3 3 6; 4 6 8; 4 7 7]
%       dbmean(X,1)
%       dbmean(X,2)

y = db2pow(x);
y = mean(y,varargin{:});
y = pow2db(y);

end