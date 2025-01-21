function L = Get_Labeling(m,type)
% L=Get_Labeling(m,type) returns the decimal representation of some
% well-known binary labelings. L is a length-M column vector and m is the
% number of bits per symbol, i.e., m=log2(M);
% The supported labelings are 'BRGC', 'NBC', and 'AGC'
%
% Code originally written by Fredrik Br�nnstr�m.
% Modified by Alex Alvarado, June 2014
%

if m==1
    L=[0 1]';
else
    M=2^m;
    switch type
        case 'BRGC'
            L=zeros(M,m);
            L(1:M/2,2:m)=Get_Labeling(m-1,type);
            L(M/2+1:M,2:m)=flipud(L(1:M/2,2:m));
            L(M/2+1:M,1)=1;
        case 'NBC'
            L=fliplr(de2bi(0:M-1));
        case 'AGC'  % This is the AGC for PAM
            L=zeros(M,m);
            L(1:M/2,2:m)=Get_Labeling(m-1,type);
            L(M/2+1:M,2:m)=flipud(L(1:M/2,2:m));
            L(2:2:M,1)=1;
            L(M/2+1:M,:)=not(L(M/2+1:M,:));
        otherwise
            error('Only ''BRGC'', ''NBC'', and ''AGC'' are supported and here type=''%s''',type);
    end
end