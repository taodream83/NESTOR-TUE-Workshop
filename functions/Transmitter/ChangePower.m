function [OutSignal,ScalingFactor,InputPower,InputPower_dBm,...
    OutputPower_Watts] = ChangePower(InSignal,DesiredPower_dBm)

% ChangePower - Vinicius Oliari - 24/01/2019
%
% This function changes the power of an input signal InSignal to a
% desired power DesiredPower_dBm in dBm. The new signal with power
% DesiredPower_dBm is the output of the function OutSignal. In addition,
% the function converts the power in dBm to Watts.
%
% %%%% Remarks %%%%
%
% It is not necessary to attribute any other output variable than the
% OutSignal (the function can be only used to normalize the signal)
%
% It is possible to use the function only for obtaining the power of the
% input signal, which can be made by calling the function as
% [~,~,InputPower,~,~] = ChangePower(InSignal,1); in Watts or
% [~,~,~,InputPower_dBm,~] = ChangePower(InSignal,1); in dBm
%
% It is possible to use the function only for convertion from dBm to
% Watts, which can be made by calling the function as
% [~,~,OutputPower_Watts] = ChangePower(1,DesiredPower_dBm);
%
% %%%%% INPUTS %%%%%
%
% InSignal: Input Signal -> signal (NxT1) that will have its power
% changed to DesiredPower_dBm in dBm. T1 is the time dimension.
%
% DesiredPower_dBm: Desired Power in dBm -> the desired power of the
% output signal (a scaled version of the input signal)
%
% %%%%% OUTPUT %%%%%
%
% OutSignal: Output Signal -> Desired output signal (NxT1), a scaled
% version of InSignal, with new % power given by DesiredPower_dBm in dBm.
%
% ScalingFactor: Scaling Factor -> Multiplication factor that will be
% applied to InSignal to change its power to DesiredPower_dBm.
%
% InputPower: Input Power -> Power of the Input Signal in Watts.
%
% InputPower_dBm: Input Power in dBm -> Power of the Input Signal in dBm.
%
% OutputPower_Watts: Output Power in Watts -> Conversion of
% DesiredPower_dBm in dBm to Watts.


% Power, in Watts, of the input signal
InputPower = mean(sum(abs(InSignal).^2,1));

% Input power in dBm
InputPower_dBm  = 10*log10(InputPower)+30;

% Calculates the Scaling Factor
ScalingFactor   = 10^((DesiredPower_dBm-InputPower_dBm)/20);

% Rescales the input signal to obtain a DesiredPower_dBm dBm power.
OutSignal       = ScalingFactor*InSignal;

% Converts DesiredPower_dBm from dBm to Watts.
OutputPower_Watts = 10^((DesiredPower_dBm-30)/10);

end

