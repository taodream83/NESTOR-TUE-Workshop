function  [Sout,Pout] = ProductDecoder(S,P)
%--------------------------------------------------------------------------
% This function performs the encoding on PC-FEC with BCH component codes.
% Jan. 2020, Initial PC decoder: Alireza Sheikh
% Jan. 2021, Modification: Yunus Can Gultekin
%--------------------------------------------------------------------------
% Input
% P: System Paramaters
% S: Signal Paramaters
%--------------------------------------------------------------------------
% Output
% P: System Paramaters
%    - .DecInfBitSeq, the decoded bit sequences
% S: Signal Paramaters - No updates
%--------------------------------------------------------------------------

Sout = S;
Pout = P;

%------ Existing parameters
s  = P.FEC.s;
ns = Pout.FEC.ns;
ks = Pout.FEC.ks;
Nblock=P.FEC.Nfecframe;
itermax=P.FEC.itermax;
%------ Decoding

% Computes TX channel indices
if mod(P.Sys.Nch,2)
    CUTidx=P.Rx.CUT.'+(P.Sys.Nch-1)/2;
else
    CUTidx1=P.Rx.CUT(P.Rx.CUT<0).'+P.Sys.Nch/2;
    CUTidx2=P.Rx.CUT(P.Rx.CUT>0).'+P.Sys.Nch/2-1;
    CUTidx=[CUTidx1 ; CUTidx2];
end

TxChIdx = (repmat(1:P.Sys.Npol,length(P.Rx.CUT),1)+CUTidx*P.Sys.Npol).'; % Transmitted channel index
TxChIdx=TxChIdx(:);


% LLR depadding to match codeword length
state_vec=P.FEC.state(TxChIdx);
DecinfBitSeq=zeros(P.Sys.Npol,ks^2);

for jj=1:P.Sys.Npol*length(P.Rx.CUT)

    state=state_vec(jj);
    LLR_vec=S.RxSym.LLR(jj,:);
    bit_vec=(sign(LLR_vec.')+1)/(2);
    bit_vec_deint=randdeintrlv(bit_vec,state);

    stack_PC_mat=reshape(bit_vec_deint,Nblock*ns,ns);

    Mat_dec_info=zeros(Nblock*ks,ks);

    for kk=1:Nblock

        vec=stack_PC_mat((kk-1)*ns+1:(kk)*ns,:);

        old_vec=vec;

        for iter=1:itermax


            for i=1:ns

                [cv,~]=step(P.FEC.Decobj,[zeros(s,1);vec(i,:)']);
                cvp=step(P.FEC.Encobj,cv);
                DHhat=sum(mod(cvp+[zeros(s,1);vec(i,:)'],2));

                if (DHhat>P.FEC.t) || isempty((find(cvp(1:s)==zeros(s,1),1)))      % bounded distance decoding
                    %                        if (DHhat>S.FEC.t) || (length(find(cvp(1:s)~=zeros(s,1))))    % bounded distance decoding

                    vec(i,:)=vec(i,:);

                else

                    vec(i,:)=cvp(s+1:end)';

                end

            end

            vec2=vec';

            for i=1:ns

                [cv,~]=step(P.FEC.Decobj,[zeros(s,1);vec2(i,:)']);
                cvp=step(P.FEC.Encobj,cv);
                DHhat=sum(mod(cvp+[zeros(s,1);vec2(i,:)'],2));

                if (DHhat>P.FEC.t) || isempty((find(cvp(1:s)==zeros(s,1),1)))      % bounded distance decoding
                    %                        if (DHhat>S.FEC.t) || (length(find(cvp(1:s)~=zeros(s,1))))    % bounded distance decoding

                    vec2(i,:)=vec2(i,:);

                else

                    vec2(i,:)=cvp(s+1:end)';

                end

            end

            vec=vec2';
            new_vec=vec;

            if isempty(find(new_vec~=old_vec,1))

                break;

            else

                old_vec=vec;

            end


        end %% End of iteration

        Mat_dec_info((kk-1)*ks+1:(kk)*ks,:)=vec(1:ks,1:ks);


    end  %% End of decoding for all frames

    DecinfBitSeq(jj,:)= Mat_dec_info(:).';

end  %% End of decoding for both polarizations

Sout.RxSym.EstimatedBitsHD = DecinfBitSeq;

if Pout.FEC.pad~=0
    Sout.CodeIntBits=Sout.CodeIntBits(:,1:end-Pout.FEC.pad);                      % Zero-padding to map to the closest number of symbols containing the coded bit sequence
end

