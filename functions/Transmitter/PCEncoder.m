function [Sout,Pout] = PCEncoder(S,P)

% This function performs the encoding on PC FEC with BCH component codes

% %%%% INPUTS %%%%%
%
% S: Input Signal Structure -> S.FEC.v: the component code is built in GF(2^v)
%                           -> S.FEC.t: error correction capability of the componenr code
%                           -> S.FEC.s: shorthening parameter for component code

% P: Structure with the link/fiber/sim parameters

% %%%% OUTPUT %%%%%
%
% S: Output Signal Structure  -> S.BitSeqNt: the number of symbols per polarization per channel
%                             -> S.FEC.PC_error: the number of inf. errors
%                             -> S.FEC.state: the seed for random interleaver for all channels and both polarizations
%                             -> S.BitSeq: encoded bits for all channels and both polarizations
%                             -> S.infBitSeq: inf. bits all channels and both polarizations
%%%%% ========================================================= %%%%%%

%     PC decoder - Alireza Sheikh, Jan. 2020

%%%% ========================================================== %%%%%%

v=P.FEC.v;
t=P.FEC.t;
s=P.FEC.s;

Pout=P;
n=2^(v)-1;                 % original component block length
k=2^(v)-v*t-1;             % original component information length
ns=2^(v)-1-s;              % Shorthened component block length
ks=2^(v)-v*t-1-s;          % Shorthened component information length


Nblock    	= ceil((P.Tx.Nsym*P.m)/((ks)^2));         % Number of PC blocks.
if Nblock<1
    Nblock=1;
end
Nbit          =Nblock*((ks)^2);                    % number of information bits in the tx sequence
Nsym          =ceil(Nbit/P.m);                  % update the code bit sequence due to the PC structure
Ncbit         =Nblock*((ns)^2);
Ncsym         =ceil(Ncbit/P.m);

enc = comm.BCHEncoder(n,k);
dec = comm.BCHDecoder(n,k,'NumCorrectedErrorsOutputPort',true, 'GeneratorPolynomialSource', 'Auto');

%%%%%%%%%%% ============================= Encoder ================================== %%%%%%%%
infBitSeq= zeros(P.Sys.Npol*P.Sys.Nch,Nbit);
BitSeq= zeros(P.Sys.Npol*P.Sys.Nch,Ncbit);
FECstate= zeros(1,P.Sys.Npol*P.Sys.Nch);

for ii=1:P.Sys.Nch

    for jj=1:P.Sys.Npol

        inf_bit=randi([0,1],Nblock*ks,ks);                       % Information bits
        infBitSeq((ii-1)*P.Sys.Npol+jj,:) = inf_bit(:).';

        encod_PC_block1=zeros(Nblock*ks,ns);     % PC Block
        encod_PC_block=zeros(Nblock*ns,ns);      % PC Block

        for kk=1:Nblock

            %%% Row Enc.
            for tt=1:ks
                en_temp1=inf_bit((kk-1)*ks+tt,:)';
                en_temp2=step(enc,[zeros(s,1);en_temp1])';
                encod_PC_block1((kk-1)*ks+tt,:)=en_temp2(s+1:end);
            end

            %%% Col Enc.
            for tt=1:ns
                en_temp1=encod_PC_block1((kk-1)*ks+1:(kk)*ks,tt);
                en_temp2=step(enc,[zeros(s,1);en_temp1]);
                encod_PC_block((kk-1)*ns+1:(kk)*ns,tt)=en_temp2(s+1:end);
            end


        end

        senddata=encod_PC_block(:);
        state=randi(100000000);
        FECstate((ii-1)*P.Sys.Npol+jj)=state;
        senddata=randintrlv(senddata,state);
        BitSeq((ii-1)*P.Sys.Npol+jj,:)=senddata.';

    end     % End of Pol

end         % End of Channel
Sout=S;
Sout.InfBitSeq=infBitSeq;               % Information bit sequence
Sout.BitSeq=BitSeq;                   % Coded bit sequence
Pout.Tx.Nbit=Nbit;        	      % Number of info bits per channel per polarization
Pout.Tx.Nsym=Nsym;        	      % Number of info syms per channel per polarization
Pout.FEC.Nblock=Nblock;                % Number of PC blocks in TX sequence
Pout.FEC.state=FECstate';             % FEC state
Pout.FEC.Encobj=enc;                  % Encoder object
Pout.FEC.Decobj=dec;                  % Decoder object
Pout.FEC.Ncbit=Ncbit;        	      % Number of coded bits per channel per polarization
Pout.FEC.Ncsym=Ncsym;        	      % Number of coded syms per channel per polarization
Pout.FEC.n=ns;
Pout.FEC.k=ks;

end