function [Sout,Pout] = ProductEncoder(S,P)
%--------------------------------------------------------------------------
% This function performs the encoding of PC-FEC with BCH component codes.
% Jan. 2020, Initial PC decoder: Alireza Sheikh
% Jan. 2021, Oct. 2023 Modification: Yunus Can Gultekin
%--------------------------------------------------------------------------

Sout = S;
Pout = P;

%------ Existing parameters
s      = P.FEC.s;
nc     = P.FEC.nc;
kc     = P.FEC.kc;
ns     = P.FEC.ns;
ks     = P.FEC.ks;
Nblock = P.FEC.Nfecframe;
Ncbit  = P.FEC.Ncbit;

%------ Derived parameters
enc    = comm.BCHEncoder(nc,kc);
dec    = comm.BCHDecoder(nc,kc,'NumCorrectedErrorsOutputPort',true, 'GeneratorPolynomialSource', 'Auto');

%------ Encoding
CodedBits = zeros(P.Sys.Npol*P.Sys.Nch,Ncbit);
FECstate  = zeros(1,P.Sys.Npol*P.Sys.Nch);

for ii=1:P.Sys.Nch % loop over channels
    for jj=1:P.Sys.Npol % loop over polarizations

        inf_bit = reshape(Sout.InfoBits((ii-1)*P.Sys.Npol+jj,:),Nblock*ks,ks);

        encod_PC_block1=zeros(Nblock*ks,ns);     % PC Block
        encod_PC_block=zeros(Nblock*ns,ns);      % PC Block

        for kk=1:Nblock

            %--- Row Encoding
            for tt=1:ks
                en_temp1=inf_bit((kk-1)*ks+tt,:)';
                en_temp2=step(enc,[zeros(s,1);en_temp1])';
                encod_PC_block1((kk-1)*ks+tt,:)=en_temp2(s+1:end);
            end

            %--- Column Encoding
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
        CodedBits((ii-1)*P.Sys.Npol+jj,:)=senddata.';
    end % of loop over polarizations
end % of loop over channels

Sout.CodedBits  = CodedBits;            % Coded bit sequence
Pout.FEC.state  = FECstate';            % FEC state
Pout.FEC.Encobj = enc;                  % Encoder object
Pout.FEC.Decobj = dec;                  % Decoder object
end