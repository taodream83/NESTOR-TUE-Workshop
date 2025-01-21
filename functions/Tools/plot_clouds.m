function plot_clouds(S,P)
% plot_clouds(S,P) plots conditional constellation clouds for the CUT.
% This function can he used to generate quick plots and see if the
% constellation clouds are where they should be.
% The colors are from the file better_colorwheel.m
%
% Author: Alex Alvarado
% Created: December 2019

better_colorwheel
CUT=P.Rx.CUT;
figure;

for chn = 1:length(P.Rx.CUT)
    for pol = 1:P.Sys.Npol
        TxChIdx = (P.Rx.CUT(chn)-1)*P.Sys.Npol+pol; % Tx index
        RxChIdx = (chn-1)*P.Sys.Npol+pol; % Rx index
        subplot(P.Sys.Npol,length(CUT),chn+(pol-1)*length(CUT));
        grid on;hold on;axis square;%cntr=cntr+1;
        vv=unique(S.TxIdx(TxChIdx,:));
        for ii=1:P.M
            txidx=find(S.TxIdx(TxChIdx,:)==vv(ii));
            plot(real(S.RxSym.EstimatedSymbols(RxChIdx,txidx)),imag(S.RxSym.EstimatedSymbols(RxChIdx,txidx)),'.')
        end
        ax = gca;ax.ColorOrderIndex = 1; % reset colors, and plot constellation points
        for ii=1:P.M
            plot(real(P.C(ii)),imag(P.C(ii)),'o','LineWidth',1.5,'MarkerFaceColor','y','MarkerSize',7)
        end
        title(['Channel #',num2str(P.Rx.CUT(chn)),', Pol #',num2str(pol)]);
    end
end

end