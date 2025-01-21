function Pout = ParseParameters(P)

Pout = P;

% Define constants
Pout.constant.c = 299792458; % light speed [m/s]
Pout.constant.h = 6.62607015e-34; % Plank constant (J/s)

%% Check for duplicate parameters
checkVariables(P,'Fibre.D','Fibre.beta2');
checkVariables(P,'Fibre.S','Fibre.beta3');
checkVariables(P,'Fibre.Or','Fibre.Ns','Sim.Fs');
checkVariables(P,'Sim.SSFM.dz','Sim.SSFM.nsteps');

%% Warnings & errors
% System parameters errors
if strcmpi(P.FEC.CodeType,'none') % uncoded
    if isfield(P.Tx, 'Nsym')
        if isfield(P.FEC,'Nfecframe') % still, Nfecframe is set
            warning('Preset parameter ''# of FEC frames'' will be discarded since transmission is uncoded.')
        else
        end
    else
        error('Number of TX symbols is not specified!')
    end
else % coded
    if isfield(P.FEC,'Nfecframe') && isfield(P.Tx, 'Nsym')
        warning('Preset parameter ''# of TX symbols'' will be overwritten w.r.t. the preset parameter ''# of FEC frames''.')
    elseif isfield(P.FEC,'Nfecframe') && ~isfield(P.Tx, 'Nsym')
    elseif ~isfield(P.FEC,'Nfecframe') && isfield(P.Tx, 'Nsym')
        warning('Parameter ''# of FEC frames'' will be computed w.r.t. the preset parameter ''# of TX symbols''.')
    else
        error('Neither # of TX symbols nor # of FEC frames is not specified!')
    end
end

if ~isfield(P.Sys,'lambda')
    error('Transmission wavelength not specified')
end
if ~isfield(P.Sys,'Npol')
    error('Number of polarisation channels not specified')
end
if ~isfield(P.Sys,'Nch')
    error('Number of TX WDM channels not specified')
end
if ~isfield(P.Sys,'Rs')
    error('Symbol rate not specified')
end
if ~isfield(P.Sys,'Chsp')
    error('WDM channel spacing not specified')
end
if ~isfield(P.Sys,'Pch')
    error('TX power per channel not specified')
end
if ~(isscalar(P.Sys.Pch) || length(P.Sys.Pch) == P.Sys.Nch)
    error('Launch power per channel needs to be scalar or a vector with length equal to the number of channels')
end

if ~(P.N==2 || P.N==4)
   error('Unsupported number of channel dimensions')
elseif (P.N == 4 && P.Sys.Npol==1) % To be extended when 1D transmission is available
    error('4D constellations and single polarisation transmission are incompatible')
end

%% Fibre link parameters errors
if isfield(P,'Fibre')     % Checks if fibre channel is used in the simulation

    if isfield(P.Fibre,'D')
        beta2=-P.Fibre.D*P.Sys.lambda.^2/(2*pi*Pout.constant.c)*1e-24;
    elseif isfield(P.Fibre,'beta2')
        beta2=P.Fibre.beta2*1e-27;
    else
        error('Neither of fibre chromatic dispersion parameters is specified');
    end

    if ~isfield(P.Fibre,'alpha')
        error('Fibre attenuation coefficient not specified');
    end

    if ~isfield(P.Fibre,'gamma')
        error('Fibre nonlinearity coefficient not specified');
    end

    % Birefringence/PMD inconsistencies
    if isfield(P.Fibre,'PMD')
        if ~isfield(P.Fibre,'Lcorr')
            warning('PMD correlation length not specified');
            warning('SOP correlation length is set to 100m');
            Pout.Fibre.Lcorr=0.1;
        end

    elseif isfield(P.Fibre,'Lcorr')
        warning('Lcorr parameter specified but PMD parameter not specified')
    end


    if ~isfield(P.Link,'nspans')
        error('Number of fibre spans not specified');
    end

    if ~isfield(P.Link,'spanlength')
        error('Fibre span length not specified');
    end



    % Errors and warnings related to propagation equation options
    if ~isfield(P.Fibre, 'PropModel')
        error('Fibre propagation model not specified');
    end

    if strcmp(P.Fibre.PropModel,'NLSE')

        if (isfield(P.Fibre,'Biref') && P.Fibre.Biref==1)
            warning('The NLSE mode of propagation does not support the emulation of the fibre birefringence. The birefringence option will be dismissed.');
        end
        if isfield(P.Fibre,'PMD')
            warning('The NLSE mode of propagation does not support PMD emulation. The PMD parameters will be dismissed.');
        end

    elseif ~strcmp(P.Fibre.PropModel,'Manakov')
        error('Fibre propagation model not correctly specified');
    end

    %% SSFM warnings
    if isfield(P.Sim.SSFM, 'Type') && isequal(P.Sim.SSFM.Type, 'uniform')
        if isfield(P.Sim.SSFM,'PhiMax')
            warning('PhiMax parameter specified with a uniform step size SSFM, simulation will be performed as uniform step size SSFM');
        end
    end

    %% Channel memory warning
    if (P.Tx.Nsym<2*abs(beta2)*P.Link.nspans*P.Link.spanlength*1e3*P.Sys.Rs*P.Sys.Nch*P.Sys.Chsp)
        warning('Channel memory longer or comparable to simulation time window');
    end

end


% Check channel under test (CUT) vector
if ~isfield(P.Rx, 'CUT')
 if mod(P.Sys.Nch,2)
    Pout.Rx.CUT=0;                            % Sets by default central channel as CUT
 else
    Pout.Rx.CUT=1;                            % Sets by default first channel on the right-hand side of the spectrum as CUT
 end

elseif ischar(P.Rx.CUT) && strcmp(P.Rx.CUT,'all')

    if mod(P.Sys.Nch,2)
        Pout.Rx.CUT = -(P.Sys.Nch-1)/2:(P.Sys.Nch-1)/2;      % normalized center frequencies
    else
        Pout.Rx.CUT = [-P.Sys.Nch/2:-1,1:P.Sys.Nch/2];
    end

end

if max(abs(Pout.Rx.CUT)) > P.Sys.Nch/2
    error('Channel Under Test (CUT) index incompatible with selected number of channels (Nch)');
end



%% Check GPU/toolbox availability
if P.Flags.GPU == 1
    try
        if gpuDeviceCount == 0
            Pout.Flags.GPU = 0; % For machines without GPU
            warning('No compatible GPU found, falling back to CPU');
        end
    catch
        Pout.Flags.GPU = 0; % For machines without the parallel computing toolbox
        warning('Parallel Computation Toolbox not found, falling back to CPU');
    end
end

end

function checkVariables(P,varargin)
%% Throw an error if any amount of input parameters exist simultaneously

numExist = zeros(1,length(varargin));
for n = 1:length(varargin)
    numExist(n) = issubfield(P,varargin{n});
end

if sum(numExist) > 1
    errorText = 'The following parameters should not exist together:';
    for n = 1:length(varargin)
        if numExist(n) ~= 0
            errorText = sprintf('%s, %s',errorText,varargin{n});
        end
    end
    error(errorText);
end

end

function r = issubfield(s,f)

% ISSUBFIELD tests for the presence of a field in a structure just like the standard
% Matlab ISFIELD function, except that you can also specify nested fields
% using a '.' in the fieldname. The nesting can be arbitrary deep.
%
% Use as
%   f = issubfield(s, 'fieldname')
% or as
%   f = issubfield(s, 'fieldname.subfieldname')
%
% This function returns true if the field is present and false if the field
% is not present.
%
% See also ISFIELD, GETSUBFIELD, SETSUBFIELD

% Copyright (C) 2005-2013, Robert Oostenveld
%
% This file is part of FieldTrip, see http://www.fieldtriptoolbox.org
% for the documentation and details.
%
%    FieldTrip is free software: you can redistribute it and/or modify
%    it under the terms of the GNU General Public License as published by
%    the Free Software Foundation, either version 3 of the License, or
%    (at your option) any later version.
%
%    FieldTrip is distributed in the hope that it will be useful,
%    but WITHOUT ANY WARRANTY; without even the implied warranty of
%    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%    GNU General Public License for more details.
%
%    You should have received a copy of the GNU General Public License
%    along with FieldTrip. If not, see <http://www.gnu.org/licenses/>.
%

if isempty(f) || isempty(s)
    r = false;
else
    t = textscan(f,'%s','delimiter','.');
    t = t{1};
    r = true;
    for k = 1:numel(t)
        if isfield(s, t{k})
            s = s.(t{k});
        else
            r = false;
            return;
        end
    end
end

end