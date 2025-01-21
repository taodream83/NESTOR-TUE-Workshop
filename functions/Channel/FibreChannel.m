function [Sout,Pout] = FibreChannel(S,P)
% Two axis Manakov
% Symmeterised Split step Manakov equation solver, includes 2nd and 3rd order
% dispersion, PMD, SPM and XPM.
%
% SignalOut=Manakov(SignalIn,P)
%
% Inputs:
%  SignalIn - input signal structure
%  P - Fibre parameters structure contains
%
%   .Sim.SSFM -                                            % Structure defining main parameters of the SSFM algorithm
%               -.Type     = 'uniform';                    % Split-Step Fourier Method flavours. Options are ['uniform', 'log', 'adaptive']
%               -.dz       =0.1;                           % Step size for 'uniform' implementation of SSFM [km]
%               -.nsteps   =100;                           % Set number of steps per span instead of step size (for 'uniform' and 'log' SSFM)
%               -.PhiMax   =0.03;                          % Maximum nonlinear phase rotation allowed per step (only for 'adaptive' SSFM) [rad]
%
%
%
%  .Link -                                                 % Structure containing main optical link parameters
%               .spanlength - fibre length (m)
%
%  .Sys                                                    % Structure containing main optical system parameters
%      .lambda  - reference wavelength (m)                 % Central lambda value
%      .Plch    - Average transmitted power per channel (dBm)
%
%  .Fibre -
%        .PropModel - defines the propagation model adopted for the simulation of fibre propagation. Options are 'Manakov' for Manakov equation and 'NLSE' for the nonlinear Schroedinger equation. Default option is 'Manakov'.
%        .alpha     - fibre attenuation (Np/m)
%        .D         - dispersion parameter at reference wavelength (s/m^2)
%        .S         - dispersion slope at reference wavelength (s/m^3)
%        .PMD       - PMD parameter at reference wavelength (s/m^0.5)
%        .gamma      - nonlinear parameter (/W/m)
%        .CorrLength - distance over which the signal polarisation state loses memory (m)
%        .Biref      - Flag determining whether or not to simulate fibre birefringence in a Monte-Carlo fashion. This option is only used for the 'NLSE' propagation model. Options are 0 ('off') and 1 ('on'). Default option is 'off'.



%   Flags.GPU - [0 1] for GPU computation
%
% Returns:
%  Sout - output signal structure
%  Pout - output parameter structure

% Author: Gabriele Liga, January 2019

%% References:
% [1] C. R. Menyuk, "Nonlinear pulse propagation in birefringent optical fibers," IEEE Journal of Quantum Electronics, vol. 23, no. 2, (1987).
% [2] D. Marcuse et. al, "Application of the Manakov-PMD equation to studies of signal propagation in optical fibers with randomly varying birefringence", IEEE Journal of Lightwave Technology, vol. 15, no. 9, (1997).
% [3] G. P. Agrawal, "Nonlinear Fiber Optics", AP, (2001).
% [4] C. Prola, "PMD emulators and signal distortion in 2.48-Gb/s IM-DD lightwave systems", IEEE Photonics Technology Letters, vol. 9, no. 6, (1997).


%% Data conversion function for computation datatype and device selection
if P.Flags.DoublePrecision
    precision = @(x)double(x);
else
    precision = @(x)single(x);
end
if P.Flags.GPU
    dataFunc = @(x)gpuArray(precision(x));
else
    dataFunc = @(x)precision(x);
end

Pout=P;
[Np,Nt] = size(S.Et);                  % Number of polarisations

%% Precalculation section
Sout = MakeTimeFrequencyArray(S);  % Frequency array [Hz]
WW=2*pi*Sout.FF;                     % Omega vector [rad/s]


a=dataFunc(P.Fibre.alpha);         % Np/m
gamma=dataFunc(P.Fibre.gamma);
spanlength=dataFunc(P.Link.spanlength);


%% Defines step size distribution
switch P.Sim.SSFM.Type

    case 'uniform'     % parameter consistency to be checked
        if isfield(P.Sim.SSFM, 'nsteps')
            Nz=dataFunc(P.Sim.SSFM.nsteps);
            dz=spanlength/Nz;       % Step size

        else
            dz=P.Sim.SSFM.dz;
            Nz=ceil(spanlength/dz);
            dz=spanlength/Nz;       % Rounded up step size
        end

    case 'log'
        Nz=P.Sim.SSFM.nsteps;
        dz=spanlength/sum(logspace(0,a*spanlength/log(10),Nz))*logspace(0,a*spanlength/log(10),Nz); % Vector of logarithmically-spaced step sizes
        Pout.Sim.SSFM.dz_vec=dz;

    case 'adaptive'
        NLphiMax=P.Sim.SSFM.NLphiMax;    % [rad]
        Pout.Sim.SSFM.dz_vec=[];
end

%    Effective NL length
if ~isequal(P.Sim.SSFM.Type,'adaptive')
    if a == 0
        dzEff = dz;
    else
        dzEff = (1-exp(-a*dz))/a;
    end
end


% Computes dispersion parameters

% First-order dispersion
if isfield(P.Fibre,'D')
    B2=-P.Fibre.D*P.Sys.lambda.^2/(2*pi*P.constant.c);
else
    B2=P.Fibre.beta2;
end

% Second-order dispersion
if isfield(P.Fibre,'S') %&& (P.Fibre.S~=0)
    B3=((P.Sys.lambda/(2*pi*P.constant.c))^2)*(-4*pi*P.constant.c*B2/P.Sys.lambda+(P.Sys.lambda^2)*P.Fibre.S);
elseif isfield(P.Fibre,'beta3')
    B3=P.Fibre.beta3;
else
    B3=0;
end

d=dataFunc(1i*(B2/2*WW.^2+B3/6*WW.^3));

%% Moves DATA on GPU if GPU computation selected
Et=dataFunc(S.Et);
Ef=ifft(Et,[],2);                                            % Fourier Transform to frequency domain
Ex=Ef(1,:);
if Np==1   % Single pol case
    Ey = zeros(1,Nt);
else       % Dual pol case
    Ey=Ef(2,:);
end

if ~isfield(P.Fibre,'ISRS')    % No ISRS considered



    %% PMD case   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if isfield(P.Fibre, 'PMD') || isfield(P.Fibre, 'SOProt')   % In the PMD case only uniform and log SSFM approach allowed so far

        if isequal(P.Sim.SSFM.Type, 'adaptive')
            error('Adaptive step size with PMD not yet implemented...');
        end

        %% Recalculate step distribution based on PMD section size
        if size(dz)==1
            dzcum=(dz:dz:spanlength);
        else
            dzcum=double(single(cumsum(dz)));
        end
        dzcum=unique([dzcum,spanlength]);

        %
        Lcorr=dataFunc(P.Fibre.Lcorr);
        PMDsec=unique([Lcorr:Lcorr:(floor(spanlength/Lcorr)*Lcorr),spanlength]);
        PMDsteps=length(PMDsec);
        dzcum=unique([dzcum,PMDsec]);
        dz=diff([0,dzcum]);           % Updated step size distribution
        %Recompute effective NL length
        if a == 0
            dzEff = dz;
        else
            dzEff = (1-exp(-a*dz))/a;
        end

        Nz=length(dz);
        Pout.Sim.SSFM.dz_vec=dz;
        Pout.Fibre.PMDsec=gather(PMDsec);

        % Rotation coefficients for PMD:
        anglesInit = randn(4,PMDsteps);
        anglesNormalized = sqrt(sum(anglesInit.^2,1));
        anglesNormalized = anglesInit./repmat(anglesNormalized,size(anglesInit,1),1);
        theta = acos(anglesNormalized(1,:));
        axis = anglesNormalized(2:4,:)./repmat(sin(theta),3,1);
        anglesUniform = axis.*repmat(theta,3,1); % these angles scatter the SOP uniformly

        Pout.Fibre.Angles=anglesUniform;

        pauliSpins(:,:,1) = [1,  0;...
                             0, -1];

        pauliSpins(:,:,2) = [0,  1;...
                             1,  0];

        pauliSpins(:,:,3) = [ 0, -1i;...
                             1i,  0];


        rotation = @(angles,pauliSpins) expm(-1i*(angles(1)*pauliSpins(:,:,1)+...
            angles(2)*pauliSpins(:,:,2)+angles(3)*pauliSpins(:,:,3)));

        %% DGD calculations
        if isfield(P.Fibre, 'PMD') && (P.Fibre.PMD~=0)                     %% PMD on
            mean_DGD = P.Fibre.PMD*sqrt(P.Link.spanlength);                    % mean DGD (ps)
            DGD_sec_mean= mean_DGD/(0.9213*sqrt(PMDsteps));                        % mean DGD per step
            DGD_sec_std = DGD_sec_mean/5;                                          % DGD std per step (see [4])
            DGD_sec = dataFunc(DGD_sec_mean+DGD_sec_std*randn(1,PMDsteps));                  % DGD per section
            if (PMDsec(end)-PMDsec(end-1))<Lcorr
                DGD_sec(end)=DGD_sec(end)*(PMDsec(end)-PMDsec(end-1))/Lcorr;
            end

        else
            DGD_sec=zeros(1,PMDsteps);
        end

        Pout.Fibre.DGD=gather(DGD_sec);
        dsec=1;

        for nn=1:Nz

            D = exp(dz(nn)*d);                           % Dispersion operator
            NLFactor = 8/9*dzEff(nn)*1i*gamma;

            % 1st order PMD operator (Random DGDs per step)
            DGD_dz=DGD_sec(dsec)*(dz(nn)/Lcorr);

            PMD =exp(-1i*(DGD_dz/2)*WW);
            Dx = PMD.*D;
            Dy=conj(PMD).*D;

            Ex=Ex.*Dx;                                                   % Apply Dispersion operator
            Ey=Ey.*Dy;                                                   % Apply Dispersion operator

            Ex=fft(Ex);                                            % Fourier Transform to time domain
            Ey=fft(Ey);                                            % Fourier Transform to time domain

            N=NLFactor*(abs(Ex).^2 + abs(Ey).^2);                     % Nonlinear operator SPM & XPM

            Ex=Ex.*exp(N-dz(nn)*a/2);                                   % Apply Nonlinear and loss operators
            Ey=Ey.*exp(N-dz(nn)*a/2);                                   % Apply Nonlinear and loss operators


            Ex=ifft(Ex);                                           % Fourier Transform to frequency domain
            Ey=ifft(Ey);                                           % Fourier Transform to frequency domain

            % Rotation on PSPs frame of reference
            if dzcum(nn)==(PMDsec(dsec))           % Flag to indicate PMD section
                J = rotation(anglesUniform(:,dsec), pauliSpins);               % Make polarisation scrambling matrix
                Ef=J*[Ex;Ey];                                                  % Apply random polarisation rotation
                Ex=Ef(1,:);
                Ey=Ef(2,:);
                dsec=dsec+1;
            end
        end


        %% No PMD and no SOP rotations %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    else

        if Np==1       %Single pol

            if isequal(P.Sim.SSFM.Type, 'adaptive')
                zcum=0;
                stop=0;
                L=P.Link.spanlength;
                PeakPow=max(sum(abs(Et).^2,1));
                dz=NLphiMax/(gamma*PeakPow);

                while ~stop

                    if (L-zcum)>dz
                        zcum=zcum+dz;
                    else
                        dz=L-zcum;
                        zcum=L;
                        stop=1;
                    end
                    Pout.Sim.SSFM.dz_vec=[Pout.Sim.SSFM.dz_vec,dz];
                    dzEff = (1-exp(-a*dz))/a;
                    NLFactor = 8/9*dzEff*1i*gamma;
                    D = exp(dz*d);                           % Dispersion operator

                    Ex=Ex.*D;                                                   % Apply Dispersion operator

                    Ex=fft(Ex);                                            % Fourier Transform to time domain
                    PeakPow=max(abs(Ex).^2);

                    N=NLFactor*(abs(Ex).^2);                      % Nonlinear operator SPM & XPM

                    Ex=Ex.*exp(N-dz*a/2);                                   % Apply Nonlinear and loss operators

                    Ex=ifft(Ex);                                           % Fourier Transform to frequency domain

                    % Compute new step size
                    dz=NLphiMax/(gamma*PeakPow);

                end

            elseif size(dz)==1    % Uniform step-size
                D = exp(dz*d);                           % Half-step dispersion operator
                NLFactor = 8/9*dzEff*1i*gamma;
                for nn=1:Nz

                    Ex=Ex.*D;                                                   % Apply Dispersion operator

                    Ex=fft(Ex);                                            % Fourier Transform to time domain

                    N=NLFactor*(abs(Ex).^2);                      % Nonlinear operator SPM & XPM

                    Ex=Ex.*exp(N-dz*a/2);                                   % Apply Nonlinear and loss operators

                    Ex=ifft(Ex);                                           % Fourier Transform to frequency domain
                end

            else    %% Log-step
                for nn=1:Nz

                    D = exp(dz(nn)*d);                                          %  Dispersion operator

                    Ex=Ex.*D;                                                   % Apply Dispersion operator

                    Ex=fft(Ex);                                                 % Fourier Transform to time domain

                    N=8/9*dzEff(nn)*1i*gamma*(abs(Ex).^2);                      % Nonlinear operator SPM & XPM

                    Ex=Ex.*exp(N-dz(nn)*a/2);                                   % Apply Nonlinear and loss operators

                    Ex=ifft(Ex);                                           % Fourier Transform to frequency domain

                end

            end

        else   %Dual pol
            if isequal(P.Sim.SSFM.Type, 'adaptive')
                zcum=0;
                stop=0;
                L=P.Link.spanlength;
                PeakPow=max(sum(abs(Et).^2,1));
                dz=NLphiMax/(gamma*PeakPow);
                while ~stop

                    if (L-zcum)>dz
                        zcum=zcum+dz;
                    else
                        dz=L-zcum;
                        zcum=L;
                        stop=1;
                    end
                    Pout.Sim.SSFM.dz_vec=[Pout.Sim.SSFM.dz_vec,dz];

                    dzEff = (1-exp(-a*dz))/a;
                    NLFactor = 8/9*dzEff*1i*gamma;
                    D = exp(dz*d);                           % Dispersion operator

                    Ex=Ex.*D;                                                   % Apply Dispersion operator
                    Ey=Ey.*D;                                                   % Apply Dispersion operator

                    Ex=fft(Ex);                                            % Fourier Transform to time domain
                    Ey=fft(Ey);                                            % Fourier Transform to time domain
                    PeakPow=max(abs(Ex).^2+abs(Ey).^2);
                    N=NLFactor*(abs(Ex).^2 + abs(Ey).^2);                      % Nonlinear operator SPM & XPM

                    Ex=Ex.*exp(N-dz*a/2);                                   % Apply Nonlinear and loss operators
                    Ey=Ey.*exp(N-dz*a/2);                                   % Apply Nonlinear and loss operators

                    Ex=ifft(Ex);                                           % Fourier Transform to frequency domain
                    Ey=ifft(Ey);                                           % Fourier Transform to frequency domain

                    % Compute new step size
                    dz=NLphiMax/(gamma*PeakPow);

                end

            elseif size(dz)==1   % Uniform step-size
                D = exp(dz*d);                           % Half-step dispersion operator
                NLFactor = 8/9*dzEff*1i*gamma;
                for nn=1:Nz

                    Ex=Ex.*D;                                                   % Apply Dispersion operator
                    Ey=Ey.*D;                                                   % Apply Dispersion operator

                    Ex=fft(Ex);                                            % Fourier Transform to time domain
                    Ey=fft(Ey);                                            % Fourier Transform to time domain

                    N=NLFactor*(abs(Ex).^2 + abs(Ey).^2);                      % Nonlinear operator SPM & XPM

                    Ex=Ex.*exp(N-dz*a/2);                                   % Apply Nonlinear and loss operators
                    Ey=Ey.*exp(N-dz*a/2);                                   % Apply Nonlinear and loss operators

                    Ex=ifft(Ex);                                           % Fourier Transform to frequency domain
                    Ey=ifft(Ey);                                           % Fourier Transform to frequency domain

                end
            else  % Log-step
                for nn=1:Nz
                    D = exp(dz(nn)*d);                                         %  Dispersion operator

                    Ex=Ex.*D;                                                 % Apply Dispersion operator
                    Ey=Ey.*D;                                                 % Apply Dispersion operator

                    Ex=fft(Ex);                                            % Fourier Transform to time domain
                    Ey=fft(Ey);                                            % Fourier Transform to time domain

                    N=8/9*dzEff(nn)*1i*gamma*(abs(Ex).^2 + abs(Ey).^2);                      % Nonlinear operator SPM & XPM

                    Ex=Ex.*exp(N-dz(nn)*a/2);                                   % Apply Nonlinear and loss operators
                    Ey=Ey.*exp(N-dz(nn)*a/2);                                   % Apply Nonlinear and loss operators

                    Ex=ifft(Ex);                                           % Fourier Transform to frequency domain
                    Ey=ifft(Ey);                                           % Fourier Transform to frequency domain

                end

            end

        end
    end

    %% ISRS case %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

else                     % ISRS case currently run only without PMD and in the adaptive step mode
    zcum=0;
    stop=0;
    L=P.Link.spanlength;
    rr=P.Fibre.Cr*FF;
    PeakPow=max(sum(abs(Et).^2,1));
    dz=NLphiMax/(gamma*PeakPow);
    Ptot=mean(sum(abs(Et).^2,1));

    while ~stop
        if (L-zcum)>dz
            zcum=zcum+dz;
        else
            dz=L-zcum;
            zcum=L;
            stop=1;
        end

        dzEff = (1-exp(-a*dz))/a;
        NLFactor = 8/9*dzEff*1i*gamma;
        RR=rr*Ptot*dzEff;
        D = exp(dz*d);                           % Dispersion operator

        Ex=Ex.*D.*exp(-dz*a/2+RR/2);                                                   % Apply Dispersion operator
        Ey=Ey.*D.*exp(-dz*a/2+RR/2);                                                   % Apply Dispersion operator

        Ex=fft(Ex);                                            % Fourier Transform to time domain
        Ey=fft(Ey);                                            % Fourier Transform to time domain
        PeakPow=max(abs(Ex).^2+abs(Ey).^2);
        N=NLFactor*(abs(Ex).^2 + abs(Ey).^2);                      % Nonlinear operator SPM & XPM

        Ex=Ex.*exp(N);                                   % Apply Nonlinear and loss operators
        Ey=Ey.*exp(N);                                   % Apply Nonlinear and loss operators

        Ptot=mean((abs(Ex).^2+abs(Ey).^2));

        Ex=ifft(Ex);                                           % Fourier Transform to frequency domain
        Ey=ifft(Ey);                                           % Fourier Transform to frequency domain

        % Compute step size
        dz=NLphiMax/(gamma*PeakPow);
        Pout.Sim.SSFM.dz_vec=[Pout.Sim.SSFM.dz_vec,dz];

    end

end


%Pout.Sim.SSFM.dz_vec=dz;

%% Gathers and returns optical field
if Np==1
    Ef = Ex;
else
    Ef = [Ex;Ey];
end
Sout.Et=fft(Ef,[],2);                              % Fourier Transform to time domain
Sout.Et = double(gather(Sout.Et));                 % GPU --> CPU if required
end


%% Input/Output Verification
% function P=IOVerify(P)
% warningnum = 0;
% if isfield(P,'Att')&&P.Fibre.alpha>1e-3, warning(['Attenuation value is large: ' num2str(P.Fibre.alpha) ' Np/m']), warningnum=warningnum+1; end
% if isfield(P,'D')&&abs(P.D)>1e-3, warning(['Dispersion value is large: ' num2str(P.D) ' s/m^2']), warningnum=warningnum+1; end
% if isfield(P,'Gamma')&&abs(P.Gamma)>1e-1, warning(['Gamma value is large: ' num2str(P.Gamma) ' W/m']), warningnum=warningnum+1; end
%
% %% TODO
% % PMD, S, dz, Length
% % if P.PMD>1e-1, warning(['PMD value is large: ' num2str(P.PMD) ' s/m^0.5']), warningnum=warningnum+1; end
%
% %% ENDTODO
%
% if (P.RefWavelength<1200e-9)||(P.RefWavelength>1700e-9), warning(['wavelength beyond fiber transmission window: ' num2str(P.RefWavelength) ' m']), warningnum=warningnum+1; end
%
% if (isfield(P,'GPU')&&P.GPU==1&&gpuDeviceCount<1)
%     warning('This machine does not support a GPU; Manakov will be computed on the CPU')
%     P.GPU=0;
%     warningnum=warningnum+1;
% end
%
% if (warningnum&&isfield(P,'verbose')&&(P.verbose>=1)), fprintf(['<strong>Manakov warnings: ' num2str(warningnum) '\n</strong>']), end
% end