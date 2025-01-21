function [Sout,Pout] = WDM_DeMux(S,P)
%% WDM COMPLEX SIGNALS DEMODULATOR
% Demultipex channels and Resamples signals to sampling rate of the matched filter
% Input:
% S     -Signal structure
% P     -Filter structure
%
% Output:
% Sout  - Output signal structure
% Pout  - Output filter structure
%
% Author: Sebastiaan Goossens, March 2019
Pout=P;
Pout.Filter.Fsin=P.Sim.Fs;
Pout.Filter.Fsout=P.Rx.Ns*P.Sys.Rs;
Sout=S;
Sout.Et=[];
%Sout.Et=zeros(P.Sys.Nch*P.Sys.Npol,length(S.Et));
if mod(P.Sys.Nch,2)
    Pout.Rx.Fchan=P.Sys.Fchan(P.Rx.CUT+(P.Sys.Nch-1)/2+1);
else
    Fchan1=P.Sys.Fchan(P.Rx.CUT(P.Rx.CUT<0)+P.Sys.Nch/2+1);
    Fchan2=P.Sys.Fchan(P.Rx.CUT(P.Rx.CUT>0)+P.Sys.Nch/2);
    Pout.Rx.Fchan=[Fchan1 Fchan2];
end

for cc = 1:length(P.Rx.CUT)
    Stemp=S;
    CompExp=exp(2j*pi*Pout.Rx.Fchan(cc)*S.TT);
    Stemp.Et=S.Et.*repmat(CompExp,P.Sys.Npol,1);
    Stemp=OpticalFilter(Stemp,P);      % Filters out channel
    Stemp=Resampling(Stemp,Pout);   % Resamples channel
    Sout.Et((cc-1)*P.Sys.Npol+1:cc*P.Sys.Npol,:)=Stemp.Et;
end

Pout.Rx.CUT=P.Rx.CUT;
Sout.Ns=P.Rx.Ns;
Sout.Fs=Pout.Filter.Fsout;
Sout.BaseBand=1;
Sout=MakeTimeFrequencyArray(Sout);

end
