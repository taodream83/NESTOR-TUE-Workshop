function Sout = DBP(S,P)
%DBP  Digital backpropagation function.
%   Sout = DBP(S,P) returns ...
%
%% Function running DBP with adaptive PMD sections
% Inputs
% - S.Et               Signal to be backpropagated
%   P.DBP.Mode:         'Ideal'- Perfect reversal of the fibre channel nonlinearity (see wiki)
%                         'Pragmatic'-   Reverses the fibre channel with a limited number of steps
%   P.DBP.PMD            Flag to enable PMD reversal
%% Parameters for ideal DBP
% Looks for fibre parameters for ideal inversion of the fibre channel
%P.Fibre.alpha                                       %Loss parameter  [Np/m]
%P.Fibre.beta2                                       %Second order dispersion coefficient [s^2/m]
%P.Fibre.beta3                                       %Third order dispersion coefficient  [s^3/m]
%P.Fibre.gamma                                       %Nonlinear parameter  [W/m]
%P.Fibre.PMD                                         %Polarization mode dispersion [s/sqrt(m)]
%P.Fibre.D                                           %Dispersion parameter at reference wavelength [s/m^2]
%P.Fibre.S
% P.Link
%        -.spanlength                                  % Span length [m]
%        -.nspans                                      % Number of spans
%        -.G                                           % Amplifier gains [dB]

%P.Sys.Pch                                          % Power per channel   [dBm]
%P.Sys.Lambda                                      % Reference wavelength [m]

%  And if P.DBP.PMD is also enabled
% - P.Fibre.Angles     Set of angles for backward PMD sections
% - P.Fibre.DGD        Set of DGD values for backward PMD sections
% - P.Fibre.
%% Parameters for pragmatic DBP
%   P.DBP.Nsteps
%   P.DBP.
% - P.Steps
%% Outputs
% - Sout.Et        Backpropagated signal
%

%% Author: Gabriele Liga
%% Created: Dec. 2019

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


%% Convert SI to Optical Units
%P=Opt2SI(1,P);
c=299792458;  % m/s

%% Precalculation section
Sf = MakeTimeFrequencyArray(S);          % Frequency array [THz]

%%
WW=2*pi*Sf.FF;

% Chromatic dispersion and third-order dispersion terms
B2=-P.Fibre.D*P.Sys.lambda.^2/(2*pi*c);   	    % s^2/m
d2=1i*B2/2*WW.^2;

if isfield(P.Fibre, 'S')
    B3=((P.Sys.lambda/(2*pi*c))^2)*(2*P.Sys.lambda*P.Fibre.D+(P.Sys.lambda^2)*P.Fibre.S);
    d3=1i*B3/6*WW.^3;
else
    d3=0;
end
d=dataFunc(d2+d3);

% Power profile
a=dataFunc(P.Fibre.alpha);                                 % Np/m
G=P.Link.G;                                                % Amplifier gain (linear units)

Nspans=dataFunc(P.Link.nspans);
Np=P.Sys.Npol;


%% Move arrays to GPU if used
Et=dataFunc(S.Et);
%dz=dataFunc(dz);
%Nz=dataFunc(Nz);
gamma=dataFunc(P.Fibre.gamma);


%% Switch optical field to frequency domain
Ef=ifft(Et,[],2);


%% Step distribution calculation from fibre parameters
if isequal(P.DBP.Mode,'Ideal')

    %% PMD or SOP rotation distributed compensation
    comp_PMD = isfield(P.Fibre, 'PMD') && isfield(P.DBP, 'PMD') && P.DBP.PMD == 1;
    comp_SOP = isfield(P.Fibre, 'SOProt') && isfield(P.DBP, 'PMD') && P.DBP.PMD == 1;
    if comp_PMD
        % PMD
        Angles=fliplr(P.Fibre.Angles);
        DGD=fliplr(P.Fibre.DGD);
        PMDsec=fliplr(P.Fibre.PMDsec);
        dsec=0;

    elseif comp_SOP
        Angles=fliplr(P.Fibre.Angles);
        PMDsec=fliplr(P.Fibre.PMDsec);
        dsec=0;
    end

    switch P.Sim.SSFM.Type
        case 'uniform'
            if comp_PMD || comp_SOP
                dz=dataFunc(fliplr(P.Sim.SSFM.dz_vec));
                Nz=dataFunc(length(dz));

            else    % No PMD compensation
                dz=P.Sim.SSFM.dz;
                Nz=dataFunc(ceil(P.Link.spanlength/dz));
                dz=dataFunc(P.Link.spanlength/Nz);       % Rounded up step size
            end

        case 'log'
            dz=dataFunc(fliplr(P.Sim.SSFM.dz_vec));
            Nz=dataFunc(length(dz));

        case 'adaptive'
            dz=dataFunc(fliplr(P.Sim.SSFM.dz_vec));
            Nz=dataFunc(length(dz));

    end

elseif isequal(P.DBP.Mode,'Pragmatic')     % Work in progress

    % Step distribution calculation from fibre parameters
    if isfield(P.DBP, 'Nsteps')           % Assumes uniform step with Nsteps/span
        Nz=dataFunc(P.DBP.Nsteps);

        if isequal(P.DBP.Type, 'uniform')
            dz=dataFunc(P.Link.spanlength/P.DBP.Nsteps);
        elseif isequal(P.DBP.Type, 'log')
            dz=dataFunc(fliplr(P.Link.spanlength/sum(logspace(0,a*P.Link.spanlength/log(10),P.DBP.Nsteps))*logspace(0,a*P.Link.spanlength/log(10),P.DBP.Nsteps)));
        else
            error('Please specify DBP step distribution (uniform or logarithmic)...');
        end


    elseif isfield(P.DBP, 'dz_vec')
        dz=P.DBP.dz_vec;
    else
        error('Select either number of steps per span or step distribution')
    end

else
    error('Undefined DBP mode. Choose between "Ideal" and "Pragmatic".');
end


if ( (Np==1) && ~(isfield(P.Fibre, 'PMD') && (P.Fibre.PMD==1)) ) % Single-pol transmission with no PMD

    if size(dz)==1
        invD=exp(-dz*d);                             % Inverse dispersion operator
        dzEff=(1-exp(-a*dz))/a;                   % Effective length of the fibre span (km)
        NLFactor = 8/9*dzEff*1i*gamma;  % NL factor

        for ss=1:Nspans
            Ef=Ef/sqrt(G);    % Renormalise before DBP

            for nn=1:Nz

                Ef=Ef.*invD;                                                % Apply Dispersion operator for the first half step

                Et=fft(Ef);                                                 % Fourier Transform to time domain

                invN=-NLFactor*(abs(Et).^2);                   % Inverse nonlinear operator

                Et=Et.*exp(invN+dz*a/2);                                     % Apply Nonlinear and loss operators at center of step

                Ef=ifft(Et);                                                % Fourier Transform to frequency domain
            end
        end
    else     % Nonuniform step size
        dzEff=(1-exp(-a*dz))/a;

        for ss=1:Nspans
            Ef=Ef/sqrt(G);    % Renormalise before DBP

            for nn=1:Nz
                invD= exp(-dz(nn)*d);                                         % Dispersion operator
                NLFactor = 8/9*dzEff(nn)*1i*P.Fibre.gamma;

                Ef=Ef.*invD;                                                  % Apply Dispersion operator

                Et=fft(Ef);                                            % Fourier Transform to time domain\
                invN=-NLFactor*(abs(Et).^2);          % Nonlinear operator SPM & XPM

                Et=Et.*exp(invN+dz(nn)*a/2);                                   % Apply Nonlinear and loss operators at center of step

                Ef=ifft(Et);                                           % Fourier Transform to frequency domain
            end

        end
    end

else  %Dual-pol TX case or single-pol but with PMD

    if comp_PMD || comp_SOP       % Undoes PMD

        rotation = @(angles,pauliSpins) expm(-1i*(angles(1)*pauliSpins(:,:,1)+...
            angles(2)*pauliSpins(:,:,2)+angles(3)*pauliSpins(:,:,3)));
        pauliSpins(:,:,1) = [1,  0;...
                             0, -1];

        pauliSpins(:,:,2) = [0,  1;...
                             1,  0];

        pauliSpins(:,:,3) = [ 0, -1i;...
                             1i,  0];


        if comp_PMD % Compensates entire PMD effect

            % Span loop
            for ss=1:Nspans
                Ef=Ef/sqrt(G);    % Renormalise before DBP
                Ex = Ef(1,:);
                Ey = Ef(2,:);

                zcum=P.Link.spanlength;
                dsec=dsec+1;

                for nn=1:Nz
                    invD= exp(-dz(nn)*d);
                    dzEff=(1-exp(-a*dz(nn)))/a;
                    NLFactor = 8/9*dzEff*1i*P.Fibre.gamma;

                    if (zcum==PMDsec(dsec))
                        Jinv=rotation(-Angles(:,dsec), pauliSpins);
                        Ef=Jinv*Ef;
                    end

                    DGDdz=DGD(dsec)*(dz(nn)/P.Fibre.Lcorr);
                    invPMD=exp(1i*(DGDdz/2)*WW);

                    Ex=invPMD.*invD.*Ef(1,:);                        % Reverse PMD vector
                    Ey=conj(invPMD).*invD.*Ef(2,:);

                    Ex=fft(Ex);                                            % Fourier Transform to time domain
                    Ey=fft(Ey);                                            % Fourier Transform to time domain

                    invN=-NLFactor*(abs(Ex).^2 + abs(Ey).^2);          % Nonlinear operator SPM & XPM
                    Ex=Ex.*exp(invN+dz(nn)*a/2);                                   % Apply Nonlinear and loss operators at center of step
                    Ey=Ey.*exp(invN+dz(nn)*a/2);                                   % Apply Nonlinear and loss operators at center of step

                    Ex=ifft(Ex);                                           % Fourier Transform to frequency domain
                    Ey=ifft(Ey);                                           % Fourier Transform to frequency domain

                    Ef=[Ex;Ey];
                    zcum=zcum-dz(nn);
                    if ~(dsec==(length(PMDsec))) && (zcum==PMDsec(dsec+1))  % Jump into next PMD section
                        dsec=dsec+1;
                    end

                end
                Ef=[Ex;Ey];
            end

        else  % Compensates only for SOP rotation

            for ss=1:Nspans
                Ef=Ef/sqrt(G);    % Renormalise before DBP
                Ex = Ef(1,:);
                Ey = Ef(2,:);

                zcum=P.Link.spanlength;
                dsec=dsec+1;

                for nn=1:Nz
                    invD= exp(-dz(nn)*d);
                    dzEff=(1-exp(-a*dz(nn)))/a;
                    NLFactor = 8/9*dzEff*1i*P.Fibre.gamma;

                    if (zcum==PMDsec(dsec))
                        Jinv=rotation(-Angles(:,dsec), pauliSpins);
                        Ef=Jinv*Ef;
                    end

                    Ex=invD.*Ef(1,:);                        % Reverse PMD vector
                    Ey=invD.*Ef(2,:);

                    Ex=fft(Ex);                                            % Fourier Transform to time domain
                    Ey=fft(Ey);                                            % Fourier Transform to time domain

                    invN=-NLFactor*(abs(Ex).^2 + abs(Ey).^2);          % Nonlinear operator SPM & XPM
                    Ex=Ex.*exp(invN+dz(nn)*a/2);                                   % Apply Nonlinear and loss operators at center of step
                    Ey=Ey.*exp(invN+dz(nn)*a/2);                                   % Apply Nonlinear and loss operators at center of step

                    Ex=ifft(Ex);                                           % Fourier Transform to frequency domain
                    Ey=ifft(Ey);                                           % Fourier Transform to frequency domain

                    Ef=[Ex;Ey];
                    zcum=zcum-dz(nn);
                    if ~(dsec==(length(PMDsec))) && (zcum==PMDsec(dsec+1))  % Jump into next PMD section
                        dsec=dsec+1;
                    end

                end
                Ef=[Ex;Ey];

            end

        end

        %% No distributed PMD compensation
    else
        if size(dz)==1
            invD=exp(-dz*d);                             % Inverse dispersion operator
            dzEff=(1-exp(-a*dz))/a;                   % Effective length of the fibre span (km)
            NLFactor = 8/9*dzEff*1i*gamma;  % NL factor

            for ss=1:Nspans
                Ef=Ef/sqrt(G);    % Renormalise before DBP
                Ex = Ef(1,:);
                Ey = Ef(2,:);


                for nn=1:Nz

                    Ex=Ex.*invD;                                                % Apply Dispersion operator for the first half step
                    Ey=Ey.*invD;                                                % Apply Dispersion operator for the first half step

                    Ex=fft(Ex);                                                 % Fourier Transform to time domain
                    Ey=fft(Ey);                                                 % Fourier Transform to time domain

                    invN=-NLFactor*(abs(Ex).^2 + abs(Ey).^2);                   % Inverse nonlinear operator

                    Ex=Ex.*exp(invN+dz*a/2);                                     % Apply Nonlinear and loss operators at center of step
                    Ey=Ey.*exp(invN+dz*a/2);                                     % Apply Nonlinear and loss operators at center of step

                    Ex=ifft(Ex);                                                % Fourier Transform to frequency domain
                    Ey=ifft(Ey);                                                % Fourier Transform to frequency domain
                end
                Ef=[Ex;Ey];
            end
        else     % Nonuniform step size
            dzEff=(1-exp(-a*dz))/a;

            for ss=1:Nspans
                Ef=Ef/sqrt(G);    % Renormalise before DBP
                Ex = Ef(1,:);
                Ey = Ef(2,:);

                for nn=1:Nz
                    invD= exp(-dz(nn)*d);                                         % Dispersion operator
                    NLFactor = 8/9*dzEff(nn)*1i*P.Fibre.gamma;

                    Ex=Ex.*invD;                                                  % Apply Dispersion operator
                    Ey=Ey.*invD;                                                   % Apply Dispersion operator

                    Ex=fft(Ex);                                            % Fourier Transform to time domain
                    Ey=fft(Ey);                                            % Fourier Transform to time domain

                    invN=-NLFactor*(abs(Ex).^2 + abs(Ey).^2);          % Nonlinear operator SPM & XPM

                    Ex=Ex.*exp(invN+dz(nn)*a/2);                                   % Apply Nonlinear and loss operators at center of step
                    Ey=Ey.*exp(invN+dz(nn)*a/2);                                   % Apply Nonlinear and loss operators at center of step

                    Ex=ifft(Ex);                                           % Fourier Transform to frequency domain
                    Ey=ifft(Ey);                                           % Fourier Transform to frequency domain
                end

                Ef=[Ex;Ey];
            end
        end

    end
end

%% Switch back to time domain
Sout=S;
Et=fft(Ef,[],2);                              % Fourier Transform to time domain
if P.Flags.GPU
    Sout.Et= double(gather(Et));
else
    Sout.Et= Et;
end

end
