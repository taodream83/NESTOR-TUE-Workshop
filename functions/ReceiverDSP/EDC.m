function [Sout] = EDC(S,P)
% EDC - Vinicius Oliari - 29/01/2019
%     - modified Gabriele Liga - 04/03/2019

% This function applies chormatic dispersion compensation to an input
% signal. All WDM chnnels are assumed to have the same propagation distance.
%
% %%%% Remarks %%%%
%
% For a dispersion managed link, it is possible to add a compensation
% factor if desired (or don't add it and put as input a different L)
%
% %%%% INPUTS %%%%%
%
% S: Input Signal Structure -> S.Et is the signal (NxT1) in which the EDC will be applied
%
% P: Structure with the link/fiber/sim parameters -> from that structure,
% it is necessary:
%       -> P.beta2: Second order dispersion coefficient beta_2
%       -> P.Link.totlength: Fiber link length to be compensated
%       -> P.Rx.EDC.Custom: Flag to enable
%       -> P.Rx.EDC.compfac (optional): compensation factor for dispersion
%       managed links
%    optional
%       -> P.Rx.CUT:  index of the channels of interest if S.BaseBand=1 (downconversion already performed by WDM_DeMux)
% %%%% OUTPUT %%%%%
%
% Sout: Output Signal Structure -> Sout.Et is the desired dispersion compensated output signal (NxT1)

c = 299792458;
Sout = S;

if isfield(S,'BaseBand') && (S.BaseBand==1) % if DEMUX has been done  (baseband EDC)
    CUT=P.Rx.CUT;
    % only the channels of interest specified by P.Rx.CUT will be individually compensated at baseband
    Sout.Et = S.Et;

    for cc = 1:length(CUT)

        % Defines frequency vector
        if isfield(S,'FF')
            WW = 2*pi*(S.FF+P.Rx.Fchan(cc));
        else
            Sout = MakeTimeFrequencyArray(S);
            WW = 2*pi*(Sout.FF+P.Rx.Fchan(cc));
        end

       

            % First-order dispersion
            if isfield(P.Fibre,'D')
                beta2 = -P.Fibre.D*P.Sys.lambda.^2/(2*pi*c);
            else
                beta2 = P.Fibre.beta2;
            end

            % Second-order dispersion
            if isfield(P.Fibre,'S')
                beta3=((P.Sys.lambda/(2*pi*c))^2)*(-4*pi*c*beta2/P.Sys.lambda+(P.Sys.lambda^2)*P.Fibre.S);
            elseif isfield(P.Fibre,'beta3')
                beta3=P.Fibre.beta3;
            else
                beta3=0;
            end

            L = P.Link.totlength;
           
            if isfield(P.Rx.EDC, 'compfac')    % Custom EDC parameters
                r = P.Rx.EDC.compfac;
            else
                r = 1;

            end

        DComp = exp(-1i*(0.5*beta2*WW.^2+beta3/6*WW.^3)*r*L);
        SIG = ifft(S.Et((cc-1)*P.Sys.Npol+1:cc*P.Sys.Npol,:),[],2); % Input signal in frequency domain
        SIG = SIG.*DComp; % Applies dispersion compensation
        Sout.Et((cc-1)*P.Sys.Npol+1:cc*P.Sys.Npol,:) = fft(SIG,[],2); % Goes back to time domain
    end

else
    % means DEMUX has not been done, thus it applies CDC jointly on the whole spectrum

    % Defines frequency vector
    if isfield(S,'FF')
        WW = 2*pi*(S.FF);
    else
        Sout = MakeTimeFrequencyArray(S);
        WW = 2*pi*(Sout.FF);
    end

    if P.Rx.EDC.Custom     % Custom EDC parameters
        beta2 = P.Rx.EDC.beta2;
        L = P.Rx.EDC.L;

        % Determines if the compensation factor was passed as input
        if isfield(P.Rx.EDC,'compfac')
            r = P.Rx.EDC.compfac;
        else
            r = 1;
        end

        % Ideal EDC
    else
        if isfield(P.Fibre,'beta2')
            beta2 = P.Fibre.beta2;
        else
            beta2 = -P.Fibre.D*P.Sys.lambda.^2/(2*pi*c);
        end

        if isfield(P.Fibre,'beta3')
        elseif isfield(P.Fibre,'S')

            beta3=((P.Sys.lambda/(2*pi*c))^2)*(2*P.Sys.lambda*P.Fibre.D+(P.Sys.lambda^2)*P.Fibre.S);
        else
            beta3=0;
        end
        L = P.Link.totlength;
        r = 1;
    end

    %FF=[FF(1:end/2),]
    DComp = exp(-1i*(0.5*beta2*WW.^2+beta3/6*WW.^3)*r*L);
    SIG = ifft(S.Et,[],2); % Input signal in frequency domain
    SIG = SIG.*DComp; % Applies the dispersion compensation
    Sout.Et = fft(SIG,[],2); % Goes back to time domain

end


end
