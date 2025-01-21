function Pout = Opt2SI(flag,P)
% Converts SI to optical units optical to be consistent with industry units.
%
%
% INPUTS:
% flag - determines the conversion way
%         1 --> Opt to SI
%        -1 --> SI to Opt
%
% P - parameter structure
%
%
% RETURNS:
% P - parameter structure
%
%
% Author: Gabriele Liga, Dec 2013

Pout = P;

Pnames = [];
if isfield(P,'Fibre');      Pnames = [Pnames;fieldnames(P.Fibre)];      end
if isfield(P,'Sim'); if isfield(P.Sim,'SSFM'); Pnames = [Pnames;fieldnames(P.Sim.SSFM)]; end; end
if isfield(P,'Sys');        Pnames = [Pnames;fieldnames(P.Sys)];        end
if isfield(P,'Link');       Pnames = [Pnames;fieldnames(P.Link)];       end
if isfield(P,'Rx');  if isfield(P.Rx,'EDC');   Pnames = [Pnames;fieldnames(P.Rx.EDC)];   end; end
N_parameters = length(Pnames);

for ii=1:N_parameters

    switch Pnames{ii}
        case 'alpha'
            if flag == -1
                Pout.Fibre.alpha = P.Fibre.alpha*(10/log(10))/1e-3;  % [Np/m] --> [dB/km]
            elseif flag == 1
                Pout.Fibre.alpha = P.Fibre.alpha/(10/log(10))/1e3;   % [dB/km] --> [Np/m]
            end

        case 'beta2'
            if isfield(P.Fibre,'beta2')
                if flag == -1
                    Pout.Fibre.beta2 = P.Fibre.beta2*((1e12)^2/1e-3);  % [s^2/m] --> [ps^2/km]
                elseif flag == 1
                    Pout.Fibre.beta2 = P.Fibre.beta2*((1e-12)^2/1e3);  % [ps^2/km] --> [s^2/m]
                end
            end
            if isfield(P.Rx.EDC,'beta2')
                if flag == -1
                    Pout.Rx.EDC.beta2 = P.Rx.EDC.beta2*((1e12)^2/1e-3);  % [s^2/m] --> [ps^2/km]
                elseif flag == 1
                    Pout.Rx.EDC.beta2 = P.Rx.EDC.beta2*((1e-12)^2/1e3);  % [ps^2/km] --> [s^2/m]
                end
            end

        case 'beta3'
            if flag == -1
                Pout.Fibre.beta3 = P.Fibre.beta3*((1e12)^3/1e-3);  % [s^3/m] --> [ps^3/km]
            elseif flag == 1
                Pout.Fibre.beta3 = P.Fibre.beta3*((1e-12)^3/1e3);  % [ps^3/km] --> [s^3/m]
            end


        case 'spanlength'
            if flag == -1
                Pout.Link.spanlength =P.Link.spanlength*1e-3;   % [m] --> [km]
            elseif flag == 1
                Pout.Link.spanlength =  P.Link.spanlength*1e3;    %  [km] --> [m]
            end

        case 'totlength'
            if flag == -1
                Pout.Link.totlength = P.Link.totlength *1e-3;   % [m] --> [km]
            elseif flag == 1
                Pout.Link.totlength  = P.Link.totlength *1e3;    %  [km] --> [m]
            end

        case 'lambda'
            if flag == -1
                Pout.Sys.lambda  = P.Sys.lambda *1e9;  % [m] --> [nm]
            elseif flag == 1
                Pout.Sys.lambda  = P.Sys.lambda *1e-9;  % [nm] --> [m]
            end

        case 'D'
            if flag == -1
                Pout.Fibre.D = P.Fibre.D*(1e12/1e9/1e-3);  % [s/m^2] --> [ps/nm/km]
            elseif flag == 1
                Pout.Fibre.D = P.Fibre.D*(1e-12/1e-9/1e3);  %  [ps/nm/km] --> [s/m^2]
            end

        case 'S'
            if flag == -1
                Pout.Fibre.S = P.Fibre.S*(1e12/(1e9)^2/1e-3);  % [s/m^3] --> [ps/nm^2/km]
            elseif flag == 1
                Pout.Fibre.S = P.Fibre.S*(1e-12/(1e-9)^2/1e3);  % [ps/nm^2/km] --> [s/m^3]
            end

        case 'PMD'
            if flag == -1
                Pout.Fibre.PMD = P.Fibre.PMD*1e12*10^(3/2);  % [s/m^0.5] --> [ps/km^0.5]
            elseif flag == 1
                Pout.Fibre.PMD = P.Fibre.PMD*1e-12*10^(-3/2);  % [ps/km^0.5]  --> [s/m^0.5]
            end

        case 'Lcorr'
            if flag == -1
                Pout.Fibre.Lcorr = P.Fibre.Lcorr*1e-3;  % [m] --> [km]
            elseif flag == 1
                Pout.Fibre.Lcorr = P.Fibre.Lcorr*1e3;  % [km]  --> [m]
            end

        case 'gamma'
            if flag == -1
                Pout.Fibre.gamma = P.Fibre.gamma/1e-3;  % [/(W*m)] --> [/(W*km)]
            elseif flag == 1
                Pout.Fibre.gamma = P.Fibre.gamma/1e3;  % [/(W*km)] --> [/(W*m)]
            end

        case 'Cr'
            if flag == -1
                Pout.Fibre.Cr = P.Fibre.Cr*1e15;   % [1/W/m/Hz] --> [1/W/km/THz]
            elseif flag == 1
                Pout.Fibre.Cr = P.Fibre.Cr*1e-15;    % [1/W/km/THz] --> [1/W/m/Hz]
            end

        case 'dz'
            if flag == -1
                Pout.Sim.SSFM.dz = P.Sim.SSFM.dz*1e-3;   % [m] --> [km]
            elseif flag == 1
                Pout.Sim.SSFM.dz = P.Sim.SSFM.dz*1e3;    % [km] --> [m]
            end

        case 'L'
            if flag == -1
                Pout.Rx.EDC.L = P.Rx.EDC.L*1e-3;   % [m] --> [km]
            elseif flag == 1
                Pout.Rx.EDC.L = P.Rx.EDC.L*1e3;    % [km] --> [m]
            end

        case 'NF'
            if flag == -1
                Pout.Link.NF = 10*log10(P.Link.NF);   % [lin] --> [dB]
            elseif flag == 1
                Pout.Link.NF = 10^(P.Link.NF/10);     % [dB] --> [lin]
            end

        case 'G'
            if flag == -1
                Pout.Link.G = 10*log10(P.Link.G);   % [lin] --> [dB]
            elseif flag == 1
                Pout.Link.G = 10^(P.Link.G/10);      % [dB] --> [lin]
            end
    end

end
end