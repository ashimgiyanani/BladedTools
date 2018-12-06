% Check the power curve of maxcap turbine

V=3:0.5:24;
[Dax,Mbeta,Mr,P,Cdax,Cp,a,theta,omr]=powercurve1('MC141',V);
plot(V,P);shg
xlabel('wind speed [m/s]');
ylabel('power [W]')
title('P-V curve LW50 - variable speed');

%%
% plot of the pitch angle
% (in a similar way the other outputs can be plotted)
figure;
plot(V,theta);
xlabel('wind speed [m/s]');
ylabel('pitch angle [degrees]')
