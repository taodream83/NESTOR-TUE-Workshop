function [Sout,Pout] = PCDecoder(S,P)

% This function performs the iterative bounded distance decoding on
% PC FEC with BCH component codes

% %%%% INPUTS %%%%%
%
% S: Input Signal Structure -> S.RxSym.LLR is the LLR vec for both polarizations
%                           -> S.FEC.state contains the seeds for random interleaver for each channel and each polarization


% P: Structure with the link/fiber/sim parameters

% %%%% OUTPUT %%%%%
%
% S: Output Signal Structure  -> S.FEC.DecinfBitSeq the decoded information bits for both polarizations
%                             -> S.FEC.PC_error the number of inf. errors

%%%%% ========================================================= %%%%%%

%     PC decoder - Alireza Sheikh, Jan. 2020

%%%% ========================================================== %%%%%%
v=P.FEC.v;
t=P.FEC.t;
s=P.FEC.s;
Nblock=P.FEC.Nblock;
itermax=P.FEC.itermax;

Pout=P;
ns=2^(v)-1-s;              % Shorthened component block length
ks=2^(v)-v*t-1-s;          % Shorthened component information length
% LLR depadding to match codeword length


state_vec=P.FEC.state((P.Sys.centre_ch_idx-1)*P.Sys.Npol+1:(P.Sys.centre_ch_idx-1)*P.Sys.Npol+2);
DecinfBitSeq=zeros(P.Sys.Npol,ks^2);

for jj=1:P.Sys.Npol

    state=state_vec(jj);
    LLR_vec=S.RxSym.LLR(P.Sys.Npol*(P.Sys.centre_ch_idx-1)+jj,:);
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

Pout.FEC.BitErr=numel(find(S.InfBitSeq((P.Sys.centre_ch_idx-1)*P.Sys.Npol+1:(P.Sys.centre_ch_idx-1)*P.Sys.Npol+2,:)~=DecinfBitSeq));
Sout.DecInfBitSeq=DecinfBitSeq;


end