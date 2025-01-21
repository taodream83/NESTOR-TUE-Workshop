function plot_constlabels(P)
% Plot constellation and labeling

h = figure;

h.Position(3) = h.Position(4); % Make sure figure fits nicely around axes
plot(P.C,'ko','MarkerFaceColor','y','MarkerSize',6)
text(real(P.C),imag(P.C),dec2bin(0:P.M-1),'FontSize',10,'HorizontalAlignment','center','VerticalAlignment','cap');

title('Constellation labels plot');
xlabel('In-Phase');
ylabel('Quadrature');

% Adjust the limits for better label readability
maxAll = max(max([real(P.C); imag(P.C)]));
limFact = 1.2;
limits = maxAll*limFact;
axis square;
axis([-limits limits -limits limits]);
end