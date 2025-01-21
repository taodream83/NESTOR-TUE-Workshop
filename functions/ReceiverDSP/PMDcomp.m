function Sout = PMDcomp(S,P)
% A function that performs PMD compensation on the input signal
%
% Input:
% S     -Signal structure
% P     -Filter structure
%
% Output:
% Sout    -Output signal structure
%
% Author: Gabriele Liga, October 2019
Sout=MakeTimeFrequencyArray(S);              % Frequency array (THz)

%% PMD COMPENSATION
rotation = @(angles,pauliSpins) expm(-1i*(angles(1)*pauliSpins(:,:,1)+...
    angles(2)*pauliSpins(:,:,2)+angles(3)*pauliSpins(:,:,3)));

% Pass PMD parameters from Manakov function to reverse pol.states evolution
Angles=fliplr(P.Fibre.Angles);
PMDsteps=size(Angles,2);

pauliSpins(:,:,1) = [1,  0;...
                     0, -1];

pauliSpins(:,:,2) = [0,  1;...
                     1,  0];

pauliSpins(:,:,3) = [ 0, -1i;...
                     1i,  0];

%% POL_EQUALISATION/PMD COMPENSATION
if isfield(P.Fibre, 'PMD') && ~(P.Fibre.PMD==0)
    % Ideal PMD Equalisation
    WW=2*pi*Sout.FF; % Hz to THz
    Ef=ifft(S.Et,[],2);

    % Pass PMD parameters from Manakov function to reverse pol.states evolution
    DGD=fliplr(P.Fibre.DGD);

    %% Reverse PMD
    for nn=1:PMDsteps
        Jinv=rotation(-Angles(:,nn), pauliSpins);
        Ef=Jinv*Ef;                                      % Derotation
        invPMD=exp(1i*(DGD(nn)/2)*WW);
        Ef(1,:)=invPMD.*Ef(1,:);                        % Reverse PMD vector
        Ef(2,:)=conj(invPMD).*Ef(2,:);
    end

    Sout.Et = fft(Ef,[],2);

    %% Compensate for SOP rotation only
elseif isfield(P.Fibre, 'SOProt')
    Et=S.Et;
    for nn=1:PMDsteps
        Jinv=rotation(-Angles(:,nn), pauliSpins);
        Et=Jinv*Et;                                          % Derotation
    end

    Sout.Et=Et;

end

end
