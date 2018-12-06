% Script to perform rainflow counting of the signals

clear; clc; close all;
% Import data from a file
StartRow = 2;
NCols = 4;
iter = (3:24)';
for idx = 1:numel(iter)
    fpath = sprintf('c:\\Users\\ashim.giyanani\\OneDrive - Windwise GmbH\\Windwise\\Projects\\Rainflow\\Fa0_V%02d00_00.csv', iter(idx))
    FormatSpec = ['%f %f %f %f'];
    % Get the header names
    fid = fopen(fpath, 'r');
    HeaderNames = textscan(fid, '%s', NCols, 'Delimiter', ';');
    fclose(fid);
    % get the data
    fid = fopen(fpath, 'r');
    Data = textscan(fid, FormatSpec, 'Delimiter', '\t', 'EmptyValue', NaN, 'HeaderLines', StartRow-1, 'returnOnError', 1);
    fclose(fid);
    % Assign data to headernames
    for k = 1:size(HeaderNames{1},1)
      eval([HeaderNames{1}{k}(1:4) '= Data{k};']);
    end

    % S-N curve parameters
    sigM = 250; % endurance limit
    m = 12; % slope of the curve
    Nk = 2e6; % number of cycles until the knee point
    T = Time(end);
    dt = mode(diff(Time));

    % Calculation of turning points and rainflow counting
    tp_Mx = sig2ext(MxHF, dt);
    tp_My = sig2ext(MyHF, dt);
    tp_Wind = sig2ext(Wind, dt);
    rf_Mx = rainflow(tp_Mx, dt);
    rf_My = rainflow(tp_My, dt);
    rf_Wind = rainflow(tp_Wind, dt);
    NCyc_Mx = rf_Mx(3,:);
    NCyc_My = rf_My(3,:);
    A_Mx = rf_Mx(1,:);
    A_My = rf_My(1,:);

    % Calculation of the damage
    Damage_Mx = sum((NCyc_Mx/Nk).*((A_Mx/sigM).^m));
    Damage_My = sum((NCyc_My/Nk).*((A_My/sigM).^m));

    % expected time to failure
    Tf_Mx = T/Damage_Mx;
    Tf_My = T/Damage_My;

    % figure, rfhist(rf_Mx,30,'ampl')
    % figure, rfhist(rf_Mx,30,'mean')
    % figure, rfmatrix(rf_Mx,30,30)

    % write the file to a comma separated file
    nc = size(rf_Mx,1) + size(rf_My,1) + size(rf_Wind,1);
    nr = max([size(rf_Mx,2), size(rf_My,2), size(rf_Wind,2)]);
    WriteData1 = nan(nr,nc);
    WriteData1(1:size(rf_Mx,2),1:5) = rf_Mx';
    WriteData1(1:size(rf_My,2),6:10) = rf_My';
    WriteData1(1:size(rf_Wind,2),11:15) = rf_Wind';
    Names = ['Mx_Amp' ';' 'Mx_Mean' ';' 'Mx_NCycles' ';' 'Mx_Tbegin' ';' 'Mx_Period', ...
                 'My_Amp' ';' 'My_Mean' ';' 'My_NCycles' ';' 'My_Tbegin' ';' 'My_Period',...
                      'Wind_Amp' ';' 'Wind_Mean' ';' 'Wind_NCycles' ';' 'Wind_Tbegin' ';' 'Wind_Period' '\r\n'];
    FormatSpec = [repmat('%4.2f;',1,14),'%4.2f \r\n'];
    fpath = sprintf('c:\\Users\\ashim.giyanani\\OneDrive - Windwise GmbH\\Windwise\\Projects\\Rainflow\\rf_Fa0_V%02d00_00.txt',iter(idx));
    fid = fopen(fpath,'w');
    fprintf(fid, Names);
    fprintf(fid, FormatSpec, WriteData1');
    fclose(fid)

    % Write the histogram data to a comma separated file
    [nCyc_Mx,bins_Mx] = rfhist(rf_Mx,20,'ampl');
    [nCyc_My,bins_My] = rfhist(rf_My,20,'ampl');
    [nCyc_Wind,bins_Wind] = rfhist(rf_Wind,20,'ampl');
    WriteData2 = [nCyc_Mx', bins_Mx', nCyc_My', bins_My', nCyc_Wind', bins_Wind'];
    Names = ['nCyc_Mx;bins_Mx;nCyc_My;bins_My;nCyc_Wind;bins_Wind\r\n'];
    FormatSpec = [repmat('%4.2f;',1,5), '%4.2f \r\n'];
    fpath = sprintf('c:\\Users\\ashim.giyanani\\OneDrive - Windwise GmbH\\Windwise\\Projects\\Rainflow\\rfhist_Fa0_V%02d00_00.txt', iter(idx));
    fid = fopen(fpath,'w');
    fprintf(fid, Names);
    fprintf(fid, FormatSpec, WriteData2');
    fclose(fid);

    % Write the rfmatrix data to a comma separated file
    [m_Mx, binx_Mx, biny_Mx] = rfmatrix(rf_Mx, 20, 20, 'ampl', 'mean');
    [m_My, binx_My, biny_My] = rfmatrix(rf_Mx, 20, 20, 'ampl', 'mean');
    [m_Wind, binx_Wind, biny_Wind] = rfmatrix(rf_Mx, 20, 20, 'ampl', 'mean');
    WriteData3 = [m_Mx, binx_Mx', biny_Mx', m_My, binx_My', biny_My', m_Wind, binx_Wind', biny_Wind'];
    Names = [repmat('matrix20x20_Mx;',1,20) 'binx_Mx;biny_Mx' repmat('matrix20x20_My;',1,20) 'binx_My;biny_My' repmat('matrix20x20_Wind;',1,20) 'binx_Wind;biny_Wind\r\n'];
    FormatSpec = [repmat('%4.2f;',1,65), '%4.2f \r\n'];
    fpath =sprintf('c:\\Users\\ashim.giyanani\\OneDrive - Windwise GmbH\\Windwise\\Projects\\Rainflow\\rfmatrix_Fa0_V%02d00_00.txt', iter(idx));
    fid = fopen(fpath, 'w');
    fprintf(fid, Names);
    fprintf(fid, FormatSpec, WriteData3');
    fclose(fid);
    clearvars WriteData1 WriteData2 WriteData3 Names FormatSpec fpath
    fprintf('[%s] Files processed %d/%d \n', datetime('now'), idx,iter(end));
end
