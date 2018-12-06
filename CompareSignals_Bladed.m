% Script to compare signals from bladed output 1 and bladed output 2 i.e. compare Bladed outputs with different settings

importBladed = 0;
bladedFile1 = 'w:\Bladed\2Benergy_MC141\combined_MC141_rev44_WG_061118_DoFactive_NoVar6p3.mat';
bladedFile2 = 'w:\Bladed\2Benergy_MC141\combined_MC141_rev44_WG_061118_DoFdeactivated_NoVar6p3.mat';
s1 = load(bladedFile1, 's');
s2 = load(bladedFile2, 's');
String = 'Blade1rotationaboutz_plane'; % Blade1x_deflection_perpendiculartorotorplane, Blade1rotationaboutz_plane
unit = 'rad'; % option of 'rad' or 'deg' or 'm'
plotSignal = 'sp';
fnames1 = fieldnames(s1.s);
fnames2 = fieldnames(s2.s);

filt_s1 = rmfield(s1.s, fnames1(find(cellfun(@isempty,strfind(fnames1, String)))));
filt_s2 = rmfield(s2.s, fnames2(find(cellfun(@isempty,strfind(fnames2, String)))));
fnames1 = fieldnames(filt_s1);
fnames2 = fieldnames(filt_s2);

for ii = 1:5:numel(fnames1)
    if strmatch(unit, 'rad')
           s1_signal = rad2deg(filt_s1.(fnames1{ii}));
           s2_signal = rad2deg(filt_s2.(fnames2{ii}));
    elseif strmatch(unit, 'deg')
           s1_signal = (filt_s1.(fnames1{ii}));
           s2_signal = (filt_s2.(fnames2{ii}));
    elseif strmatch(unit, 'm')
           s1_signal = (filt_s1.(fnames1{ii}));
           s2_signal = (filt_s2.(fnames2{ii}));
    end
    figure(ii);
    plot(s1.s.Timefromstartofsimulation, s1_signal , 'r-'); hold on;
%    plot(s2.s.Timefromstartofsimulation, s2_signal, 'b-'); hold off;
    xlabel('Time');
    ylabel(string(fnames1{ii}));
    title('Comparsion of Bladed with active doFs and deactive doFs');
%    legend('doF active', 'doF deactivated');
    legend('doF active');
        
    mean1_signal(ii) = nanmean(s1_signal); 
    mean2_signal(ii) = nanmean(s2_signal);

    max1_signal(ii) = nanmax(s1_signal); 
    max2_signal(ii) = nanmax(s2_signal);

    min1_signal(ii) = nanmin(s1_signal); 
    min2_signal(ii) = nanmin(s2_signal);
    
    % rotation when the rotor position is in azimuth at zero deg
    [idx1_azimuthpos_0deg, val1_azimuthpos_0deg] = find(round(rad2deg(s1.s.Rotorazimuthangle)) == 180);
    [idx2_azimuthpos_0deg, val1_azimuthpos_0deg] = find(round(rad2deg(s2.s.Rotorazimuthangle)) == 180);
    sp1_signal(ii) = nanmean(s1_signal(idx1_azimuthpos_0deg)); 
    sp2_signal(ii) = nanmean(s2_signal(idx2_azimuthpos_0deg));
    
    temppos{ii} = regexp(fnames1{ii}, '\d*', 'match');
    radpos(ii) = cellfun(@str2double, temppos{ii}(2));
end

figure;
if strmatch(plotSignal, 'mean')
     plotSignal1 = mean1_signal;
     plotSignal2 = mean2_signal;
elseif strmatch(plotSignal, 'min')
     plotSignal1 = min1_signal;
     plotSignal2 = min2_signal;
elseif strmatch(plotSignal, 'max')
     plotSignal1 = max1_signal;
     plotSignal2 = max2_signal;
elseif strmatch(plotSignal, 'sp')
     plotSignal1 = sp1_signal;
     plotSignal2 = sp2_signal;
end
plot(radpos, plotSignal1, 'r-'); hold on;
plot(radpos, plotSignal2, 'b-'); hold off;
xlabel('Radial position along blade (m)');
% ylabel('Rotation i.e. Twist (deg)');
ylabel('tip deflection in x-direction (flapwise) [m]')
legend('doF active', 'doF deactivated');
title('tip deflection at rotor pos = 180°');
