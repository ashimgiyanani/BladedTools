% Script to derive stats related to the step study between Flex and Bladed

% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: OTHER_FUNCTION_NAME1,  OTHER_FUNCTION_NAME2

% Author: Ashim Giyanani
% Windwise GmbH, Münster
% email address: ashimgiyanani@gmail.com
% Website: n/a
% July 2018; Last revision: 12.07.2018 18:09:06

% Steps to be followed
% Import Bladed and Flex output to be compared
% plot the plots of pitch angle, wind speed, electrical power, rotor speed, thrust force, thrust force - not tilted, tower bottom bending moment, rotor torque,

%------------- BEGIN CODE --------------

%% Files for ECD
% bladedFile1 = 'c:\Users\ashim.giyanani\OneDrive - Windwise GmbH\Windwise\Bladed\2BEnergy_MC141\combined_MC141_rev46_DoFact_Inflow_min5deg_ecd_p.mat';
% bladedFile2 = 'c:\Users\ashim.giyanani\OneDrive - Windwise GmbH\Windwise\Bladed\2BEnergy_MC141\combined_MC141_rev46_DoFact_Inflow_min5deg_ecd_n.mat';
% FlexFile1 = 'c:\Users\ashim.giyanani\OneDrive - Windwise GmbH\Windwise\Flex\ECD_rev03_126mHH\Inflow_min5deg_ed2_ecd_p.txt';
% FlexFile2 = 'c:\Users\ashim.giyanani\OneDrive - Windwise GmbH\Windwise\Flex\ECD_rev03_126mHH\Inflow_min5deg_ed2_ecd_n.txt';


clear; clc; close all;
addpath('c:\Users\ashim.giyanani\OneDrive - Windwise GmbH\Windwise\Matlab\');
addpath('C:\Users\ashim.giyanani\OneDrive - Windwise GmbH\Windwise\PostProcessing\');
addpath('c:\Users\ashim.giyanani\OneDrive - Windwise GmbH\Windwise\Matlab\src\');

PlotFig = 1;
SaveFig = 1;
Validation = 0; 
impBladed = 1; % import bladed into mat file and the workspace 0 - no import, 1 - full import
oldFile = 0; % old imports of the Bladed format were not saved to the struct variable  's'
SpecFig = 0;
FlexFlap = 0; % extract flapwise deflections from Flex Docs '1'-extract, '0' - Don't extract
FlexTwist = 1; % calculate the twist due to bend-twist coupling '1'-calculate, '0' - Don't calculate
inpWndFile = 0; % if extra files need to be used for information about Bladed
stepTime=10; 
stepStart = 55;

% Prepare the Bladed files
if impBladed == 1
  Append = 1;
  AssignVar = 1;
  Format= 'binary';
  Version = 4.6;
  prjFile = 'w:\Bladed\2Benergy_MC141\MC141_rev49_DoFact_Modal_b9tw5_ed2_ecd_p_a.$PJ';
  [matFile, ~, s] = importBladed(prjFile, Append, AssignVar, Format, Version);
else
  matFile = 'w:\Bladed\2Benergy_MC141\combined_MC141_rev49_DoFact_Modal_b9tw5_ed2_ecd_p_a.mat';
  load(matFile);
end

if inpWndFile == 1
  wndFile = 'c:\Users\ashim.giyanani\OneDrive - Windwise GmbH\Windwise\Bladed\2BEnergy_MC141\Bladed_wnd_ed2_ews_vp_a.txt';
  tq = s.Timefromstartofoutput;
  [wnd, wndMat] = FnImportWnd(wndFile, tq);
end


% Input
Delimiter = '\t'; % delimiter in the file to be read
nHead = 14; % number of lines before header names
fpath = 'c:\Users\ashim.giyanani\OneDrive - Windwise GmbH\Windwise\Flex\ECD_rev04_126mHH\Inflow_min5deg_ed2_ecd_p.txt';
String = 'VuHc|Teta2|OmRot|P_Gen|FzHR|FzYB|MyTB|MzHR|FzH2|MyYB|MyHR|MxHF|TwrClr|Mx*|My*|Mz*|Fx*|Mfric*|\wH1$|M*\w11$|Wdr*|Psi|FTip*|U*';
% String = '*Tip*|Uf*|Uk*|TwrClr|FTip*';

% Read Flex data
fid = fopen(fpath);
for i = 1:nHead+1
  tLines = fgetl(fid);                                                                            % read the header lines
end
frewind(fid);                                                                                   % return cursor one row up
NCols = sum(isspace(tLines));                                                       % number of columns in the file
DefVarNames = textscan(fid, '%s', NCols, 'delimiter', Delimiter, 'HeaderLines',14, 'MultipleDelimsAsOne', 1);                    % all the variables in the csv file
fclose(fid);

% find the indexes of the colum to be imported
Var = cellstr([DefVarNames{:}]);
ind = regexp(Var, String, 'match');                                                                % regexp to match the variable names to the string input
indf = find(not(cellfun('isempty', ind)));                                                         % find matching instances
index = [1;indf];                                                                                  % indexes to be importedd
VarNames = Var(index);                                                                             % variables that will be imported
indfs = indf-1;                                                                                    %
fmt1 = {'%f'};                                                                                     %
fmt2 = repmat({'%*f'},1, NCols-1);                                                                 %
fmt2(indfs) = {'%f'};                                                                              %
FormatSpec1 = [fmt1{:} fmt2{:} '%[^\r\n]'];                                                        % format specification for the import variables
% clearvars ind indf index indfs fmt1 fmt2 FileId NCols sub fls i j

% Extract the data into individual sensors
[fid, message] = fopen(fpath,'r');
NewData = textscan(fid, FormatSpec1, 'Delimiter', Delimiter, 'EmptyValue', NaN,...
           'headerlines', nHead+4, 'returnonError',0, 'EndOfLine', '\r\n');                 % import data
fclose(fid);
%ts = datenum([NewData{:,1}]);                                                                  % time stamp for wind speed
Data = cell2mat(NewData(1:end-1));                                                        % colelcting in one array for writing purposes
clearvars NewData
nVar = length(VarNames);
for k = 1:nVar
  eval([char(VarNames(k)) '= Data(:,k);']);                                                     % using eval for assigning variable names to the data, be careful
end

%% Extract data from Bladed files
addpath('w:\Bladed\2Benergy_MC141\')
if exist('matFile', 'var')
   filename = matFile;
else
    filename = '';
end
tempSplit = regexp(filename, '_', 'split');
prjName = regexp(tempSplit{end}, '.mat', 'split');
prjName = string(prjName(1));
if oldFile == 1
   s = load(filename)
else
    load(filename)%,...
end

%%  Calculate the bend twist and flapwise deflection parameters parameters for comparison
if FlexTwist == 1
    currentMx = [];
    currentFz = [];
    ExcelFilename = 'c:\Users\ashim.giyanani\OneDrive - Windwise GmbH\Windwise\MyProjects\Excel_Calculations\LZ69-3.0-Decrypt-Excl.Ice-Rev10.11-GH4.3-20170726_CBe.xlsx';
    queryPsi = 249.65;
    queryRotorazimuthangle = 89.71;
    [cum_locTwist, bladed_meanSig, bladed_Psi, bladed_Tpsi, flex_meanSig, flex_Psi, flex_Tpsi, flex_r, bladed_r] = FnBendTwist(ExcelFilename, matFile, Time, Psi, queryPsi, queryRotorazimuthangle, stepStart, currentMx, currentFz);
end

% read the Flex Docs to extract flapwise deflection at certain locations (actually achieved using Ftip1)
if FlexFlap == 1
    [filepath,name,ext] = fileparts(fpath);
    FlapFileName = fullfile(filepath,'Docs', [name, '.doc']);
    nHead = 561;
    % Read Flex data
    fid = fopen(FlapFileName);
    for i = 1:nHead
      tLines = fgetl(fid);                                                                            % read the header lines
    end
    NCols = 7;
    FormatSpec = ['%.3f%*f%*f%.3f%*f%*f%*f', '%[^\r\n]'];
    NewData  = textscan(fid, FormatSpec, 41, 'delimiter', Delimiter,'MultipleDelimsAsOne', 1,...
                             'EmptyValue', NaN, 'returnonError',0, 'EndOfLine', '\r\n');                    % all the variables in the csv file
    fclose(fid);
    r_Flex = cell2mat(NewData{1,1}); % radial positions in Flex Docs
    Flap_dx = cell2mat(NewData{1,2}); % flapwise deflection used in the simulations
    
end
%                   'Hublongitudinalwindspeed',...
%                  'Rotorspeed',...
%                  'Electricalpower',...
%                  'Blade2pitchangle',...
%                  'RotatinghubFx',...
%                  'RotatinghubMx',...
%                  'YawbearingFx',...
%                  'Timefromstartofoutput',...
%                  'StationaryhubFx',...
%                  'Bladeroot2Fx',...
%                  'FoundationMy',...% 'FoundationMy' 'MZT_Mbr56End1'
%                  'YawbearingMy',...
%                  'StationaryhubMy',...
%                  'YawbearingMz',...
%                  'Bladeroot1Mx',...
%                  'Bladeroot1My',...
%                  'Bladeroot1Mz',...
%                  'Bladeroot1Fx',...
%                  'Bladeroot1Fy',...
%                  'Bladeroot1Fz',...
%                  'RotatinghubMx',...
%                  'RotatinghubMy',...); 
%                  'RotatinghubMz',...);
%                  'StationaryhubMx',...);
%                  'StationaryhubMy',...);
%                  's.Winddirectionathub'
%                  s.Blade1x_deflection_perpendiculartorotorplane_69p,...  
%                  s.Blade1x_deflection_perpendiculartorotorplane_46p,...  
%                  s.Blade1rotationaboutz_plane_69p),...           
%                  s.Blade1rotationaboutz_plane_46p),...           
%                  s.Rotorazimuthangle),...                        

% Plot the results
if PlotFig == 1
    PlotVar_Flex = [...
                    VuHc,...             %1 wind speed
                    OmRot,...             %2 rotor speed
                    P_Gen,...             %3 electrical power generated
                    Teta2,...             %4 pitch angle
                    FzHR,...              %5 Thrust
                    MzHR,...              %6 Torque
                    MyTB,...              %7 tower base tilt moment
                    FzYB,...              %8 Yaw thrust
                    FzH1,...              %9 blade oop bending force
                    MyYB,...              %10 yaw bearing tilt moment
                    MyHR,...              %11 rotating hub bending moment
                    MxHF,...              %12 fixed hub yaw moment
                    MzH1,...              %13 ip bending moment
                    MyH1,...              %14 oop bending moment
                    MxH1,...              %15 torsional bending moment
                    FzH1,...              %16 in plane bending force
                    FyH1,...              %17 out of plane bending force
                    FxH1,...              %18 torsional force
                    MxHR,...              %19 bending moment
                    MzHR,...              %20 torque
                    MyHF,...              %21 fixed hub tilt moment
                    MxYB,...              %22 yaw bearing yaw moment
                    MzYB,...              %23 roll
                    MxTB,...              %24 foundations yaw tower base
                    MzTB,...              %25 foundations roll tower base
                    WdrhHc,...            %26 wind direction at hub
                    FTip1,...             % 27 flapwise deflection at tip
                    cum_locTwist(41,:)',...      %28 twist at the tip
                    cum_locTwist(28,:)',...      %29 twist time series at 0.75R
                    mod((Psi+180),360)...                     %30 Rotor position, 0 deg -> blade 1 downwards
%                     U11,...                                   % wind speed at blade down position
%                     U31,...                                   % wind speed at blade up position
                    ];
%                     U21,...                                   % wnd.spd., u-comp   R =70.57 m (left)
%                     U41,...                                   % wnd.spd., u-comp   R =70.57 m (right) 

    PlotVar_Bladed = [
                      s.Hubwindspeedmagnitude,...                              %1 wind speed                                  -> matches
                      s.Rotorspeed.*30./pi,...                                 %2 rotor speed                                 -> matches        
                      s.Electricalpower./10^3,...                              %3 electrical power generated                  -> matches  
                      rad2deg(s.Blade2pitchangle),...                          %4 pitch angle                                 -> matches       
                      s.StationaryhubFx./10^3,...                              %5 thrust                                      -> matches      
                      s.StationaryhubMx./10^3,...                              %6 torque                                      -> matches 
                      s.FoundationMy./10^3,...                                 %7 tower base tilt moment                      -> matches (NCoMmin1p6m)
                      s.YawbearingFx./10^3,...                                 %8  yaw thrust                                 -> matches
                      s.Bladeroot1Fx./10^3,...                                 %9 blade oop bending force                     -> matches
                      s.YawbearingMy./10^3,...                                 %10 yaw bearing tilt moment                    -> ### differs a bit, offset 50-100 kNm at below-rated 
                      -s.RotatinghubMy./10^3,...                                %11 rotating hub bending moment                -> matches
                      -s.StationaryhubMz./10^3,...                             %12 fixed hub yaw moment                       -> ### differs a bit offset 100kNm
                      s.Bladeroot1Mx./10^3,...                                 %13 in plane bending moment                    -> matches 
                      -s.Bladeroot1My./10^3,...                                %14 oop bending moment                         -> matches 
                      s.Bladeroot1Mz./10^3,...                                 %15  torsional moment                          -> ### differs significantly
                      s.Bladeroot1Fx./10^3,...                                 %16 in plane bending force                     -> matches  
                      -s.Bladeroot1Fy./10^3,...                                %17 out of plane bending force                 -> matches    
                      s.Bladeroot1Fz./10^3,...                                 %18 torsional force                            -> matches   
                      s.RotatinghubMz./10^3,...                                %19 bending moment                             -> matches        
                      s.RotatinghubMx./10^3,...                                %20 torque                                     -> matches               
                      s.StationaryhubMy./10^3,...                              %21 fixed hub tilt moment                       -> ### big offset around 100 kNm
                      -s.YawbearingMz./10^3,...                                %22 yaw bearing yaw moment                     -> ### differs a bit     
                      s.YawbearingMx./10^3,...                                 %23 roll                                       -> matches                   
                      -s.FoundationMz./1000,...                                %24 yaw tower base                             -> matches             
                      s.FoundationMx./1000,...                                 %25 roll tower base                            -> ### differs a bit             
                      rad2deg(s.Winddirectionathub),...                        %26 wind directin at hub                       -> input, should match                              
                      s.Blade1x_deflection_perpendiculartorotorplane_69p383,...   %27 flapwise deflection at tip
                      rad2deg(s.Blade1rotationaboutz_plane_69p383),...            %28 Torsional twist at tip
                      rad2deg(s.Blade1rotationaboutz_plane_46p420),...            %29 Torsional twist at 0.75R
                      rad2deg(s.Rotorazimuthangle),...                         %30 Rotor azimuth position
%                       wnd.wind_x0_y56,...                                      % wind speed at tip, lowest z position
%                       wnd.wind_x0_y197,...                                     % wind speed at tip, highest z position
                      ];

    plotTitle = {...
                 'u_comp wind speed, hub center',...                     %1 wind speed                  
                 'rotor speed (rpm)',...                                 %2 rotor speed                 
                 'generator power incl. mech+el. losses',...             %3 electrical power generated  
                 'pitch angle blade 2',...                               %4 pitch angle                 
                 'thrust force',...                                      %5 Thrust                      
                 'Torque',...                                            %6 Torque                      
                 'Tower base tilt moment',...                            %7 tower base tilt moment      
                 'Thrust in yaw bearing',...                             %8 Yaw thrust                  
                 'blade oop force',...                                   %9 blade oop bending force     
                 'yaw bearing tilt moment',...                           %10 yaw bearing tilt moment    
                 'RotatinghubMy',...                                     %11 rotating hub bending moment
                 'YawbearingMz',...                                      %12 fixed hub yaw moment       
                 'Bladeroot1Mx',...                                      %13 ip bending moment          
                 'Bladeroot1My',...                                      %14 oop bending moment         
                 'Bladeroot1Mz',...                                      %15 torsional bending moment   
                 'Bladeroot1Fx',...                                      %16 in plane bending force     
                 'Bladeroot1Fy',...                                      %17 out of plane bending force 
                 'Bladeroot1Fz',...                                      %18 torsional force            
                 'RotatinghubMz',...                                     %19 bending moment             
                 'RotatinghubMx',...                                     %20 torque                     
                 'StationaryhubMy',...                                   %21 fixed hub tilt moment      
                 'YawbearingMz',...                                      %22 yaw bearing yaw moment     
                 'YawbearingMx',...                                      %23 roll                       
                 'FoundationMz',...                                      %24 foundations yaw tower base 
                 'FoundationMx',...                                      %25 foundations roll tower base
                 'Wind direction',...                                    %26 wind direction at hub
                 'Flapwise deflection (r=R)',...                         %27 flapwise deflection at tip 
                 'Torsional twist (r=R)',...                             %29 Torsional twist at tip
                 'Torsional twist (r=0.75R)',...                         %30 Torsional twist at 2/3R
                 'Rotor azimuth position',...                            %31 Rotor azimuth position
%                  'wind at z = H-R',...
%                  'wind at z = H+R',...
                 };

    plotYlabels = {...
                   'hub height wind speed [m/s]',...
                   'rotor speed [rpm]',...
                   'generator power [kW]',...
                   'pitch angle [deg]',...
                   'thrust force [kN]',...
                   'Torque [kNm]',...
                   'Tower base tilt moment [kNm]',...
                   'Thrust, yaw bearing [kN]', ...
                   'blade oop fore [kN]',...
                   'YawbearingMy [kNm]',...
                   'StationaryhubMy [kNm]',...
                   'YawbearingMz [kNm]',...
                   'Bladeroot1Mx [kNm]',...
                   'Bladeroot1My [kNm]',...
                   'Bladeroot1Mz [kNm]',...
                   'Bladeroot1Fx',...
                   'Bladeroot1Fy',...
                   'Bladeroot1Fz',...
                   'RotatinghubMz [kNm]',...
                   'RotatinghubMx [kNm]',...
                   'StationaryhubMy [kNm]',...
                   'YawbearingMz [kNm]',...
                   'YawbearingMx [kNm]',...
                   'FoundationMz [kNm]',...
                   'FoundationMx [kNm]',...
                   'Wind direction [deg]',...
                   'Flapwise deflection (r=R) [m]',...
                   'Torsional twist (r=R) [deg]',...
                   'Torsional twist (r=0.75R) [deg]',...
                   'Rotor azimuth poisition [deg]',...
%                    'wind speed at z=H-R [m/s]',...
%                    'wind speed at z=H+R [m/s]',...
                   };


    fileNames = {...
                 'windspeed',...
                 'rotorspeed',...
                 'genPower',...
                 'pitch',...
                 'thrust',...
                 'torque',...
                 'towerTiltM',...
                 'thrustYawB',...
                 'Foop',...
                 'YawbearingMy',...
                 'StationaryhubMy',...
                 'YawbearingMz',...
                 'Bladeroot1Mx',...
                 'Bladeroot1My',...
                 'Bladeroot1Mz',...
                 'Bladeroot1Fx',...
                 'Bladeroot1Fy',...
                 'Bladeroot1Fz',...
                 'RotatinghubMz',...
                 'RotatinghubMx',...
                 'StationaryhubMy',...
                 'YawbearingMz',...
                 'YawbearingMx',...
                 'FoundationMz',...
                 'FoundationMx',...
                 'WindDirection',...
                 'Ftip1',...
                 'BTwist_tip1',...
                 'BTwist_0p75R',...
                 'Psi',...
%                  'windspeed_Bladebottom',...
%                  'windspeed_Bladetop',...
                 };
                 
    dt = mode(diff(Time));
    Fs = 1/dt;
    for ip = 1:size(PlotVar_Flex,2)
        p1 = figure();
        plot(s.Timefromstartofoutput, PlotVar_Bladed(:,ip), 'r.-',Time,  PlotVar_Flex(:,ip), 'k.-'); hold on;
        xx = [stepStart stepStart+stepTime];
        tmp_yy = ylim();
        yy = [tmp_yy(1) tmp_yy(2)];
        r1 = rectangle('Position', [stepStart, yy(1), stepTime, abs(yy(2)-yy(1))], 'FaceColor',[0 0 1 0.2]); hold off;
        title(plotTitle{ip});
        xlabel('Time (s)');
        ylabel(plotYlabels{ip});
        legend('Bladed','Flex', 'Location', 'southwest');
        uistack(r1,'bottom');
        grid minor
        if SaveFig == 1
           saveFolder = 'c:\Users\ashim.giyanani\OneDrive - Windwise GmbH\Windwise\Reports\Pics\';
           saveFile = sprintf('%sts_%s_%s',saveFolder, prjName, fileNames{ip});
%            cleanfigure;
           saveas(p1, saveFile,  'epsc');
%          export_fig(p1, saveFile,'-transparent','-q10', '-depsc');
%            cleanfigure;
%             matlab2tikz('interpretTickLabelsAsTex',true,'maxChunkLength',4000,'floatFormat', '%.4g', 'checkForUpdates', true,...
%                                                                                                 [saveFile,'.tex'], 'showInfo', false);
            fprintf('[%s]:Saved as (ts_%s_%s.eps) \n',datetime(),prjName,fileNames{ip})
        end % savefig       
        
%         p2 = figure();
%         [Pf, ff] = periodogram(PlotVar_Flex(:,ip), [], length(PlotVar_Flex(:,ip)), Fs, 'onesided', 'psd');
%         [Pb, fb] = periodogram(PlotVar_Bladed(:,ip), [], length(PlotVar_Bladed(:,ip)), Fs, 'onesided', 'psd');
%         loglog(ff, Pf, 'k.-', fb, Pb,'r.-');
%         legend('Flex', 'Bladed');
%         ylabel(sprintf('%s Spectra', plotYlabels{ip}),'FontSize', 12);
%         xlabel('Frequency Hz', 'FontSize', 12);
%         title(sprintf('%s Spectrum comparison', plotTitle{ip}));
%         grid on; 
%         if SaveFig == 1
%           saveas(p2, sprintf('spect_%s_%s',prjName, fileNames{ip}), 'epsc')
%         end % savefig
    end % number of plots
end % plotfig
fprintf('[%s]:in -> %s \n',datetime(), saveFolder)

p3 = figure();
plot(bladed_r, bladed_meanSig, 'r.-'); hold on;
plot(flex_r, flex_meanSig, 'k.-'); hold on;
title('Comparison of mean Twist');
xlabel('Radial Positions (m)');
ylabel('mean Twist (deg)');
legend('Bladed','Flex', 'Location', 'southwest');
grid minor; hold off;
saveFile = sprintf('%sts_%s_%s', saveFolder, prjName, 'meanTwist'); 
if SaveFig == 1
  saveas(p3, saveFile, 'epsc');
end


if SpecFig == 1
    p3 = figure();
%     plot(s.Timefromstartofoutput, s.Blade1Mz_Rootaxes_0p./1000, 'r.-',Time,  Mx11, 'k.-');
    title('Torsional moment at root');
    xlabel('Time (s)');
    ylabel('Torsional moment [kNm]');
    legend('Bladed','Flex');
    grid minor

    p4 = figure();
    plot(s.Timefromstartofoutput, s.Bladetiptotowersurface, 'r.-',Time,  TwrClr, 'k.-');
    title('Blade tip to tower distance');
    xlabel('Time (s)');
    ylabel('tip to tower deflection [m]');
    legend('Bladed','Flex');
    grid minor
   
    figure; 
    plot([s.Blade1Mz_Principalaxes_0p, s.Blade1Mz_Rootaxes_0p,s.Blade1Mz_Aerodynamicaxes_0p, s.Bladeroot1Mz]./1000)
    legend('principle', 'root_0', 'aero', 'root')
    title('Torsional moment at root');
    xlabel('Time (s)');
    ylabel('Torsional moment [kNm]');
end

if Validation == 1
    % Validation of certain non-fulfilling signals, for e.g. tower tilt moments

    % Validation of Bladed output with theoretical value
    FxHt = s.YawbearingFx; % yaw bearing thrust force
    MyHt = s.YawbearingMy; % yaw bearing tilt moment
    H = 123.9;           % hub height
    z = 3;               % height of foundation
    MyZt =  FxHt.*(H-z) + MyHt; %  tower tilt moment at  height z
    figure;
    plot(s.Timefromstartofoutput, MyZt./10^6, 'k.-', s.Timefromstartofoutput, s.FoundationMy./10^6, 'r:' );
    hold on;

    % Validation of Flex output with theoretical value
    FxHt = FzYB; % yaw bearing thrust force
    MyHt = MyYB; % yaw bearing tilt moment
    MyZt =  FxHt.*(H-z) + MyHt; %  tower tilt moment at  height z
    plot(Time, MyTB./10^3, 'b.-');
    legend('theory','Bladed', 'Flex');
    
    % Validation of the stationary hub My, Blade root My, yaw bearing My and tower base My
    % bladed signals -> s.StationaryhubMy, s.Bladeroot1My, s.YawbearingMy, s.FoundationMy
    % Flex signals -> MyHF, MyB1, MyHt, MyZt
    figure;
    plot(Time, MyHF, 'k.-', s.Timefromstartofoutput, s.StationaryhubMy./10^3, 'r.-')
    legend('Flex', 'Bladed');
    title('Hub tilt moment');

    % Calculation of yaw bearing based on forces at the hub    
    tilt = 5;  % tilt angle of the shaft wrt to horizontal [deg]
    xr_ov = 4.53; % overhang [m]
    xr_hub = -0.04;  % centre of mass of the hub [m]
    xr = (xr_ov + xr_hub); % overhang i.e. distance between the rotor centre of mass and the tower axis [m]
    zr = 2.195 - xr*sind(tilt);    % hub vertical offset from tower top [m]
    xn = -1.74; %  % centre of mass of Nacelle in the longitudinal direction wrt to tower centre axis [m]
    zn = 1.8; %   centre of mass of Nacelle in the vertical direction wrt tower top [m]
    Myr = s.StationaryhubMy; % tilting moment at the rotor [Nm]
    Mxr = s.StationaryhubMx; % Torque [Nm]
    Fyr = s.Bladeroot1Fx; % in-plane bending force on the blade [N]
    Fxr = s.StationaryhubFx; % Thrust force on the rotor
%     Fzr = s.StationaryhubFz; % when taken directly from Bladed, is in global co-ordinate so no need to convert using cosd(tilt)
    g = 9.81; % acceleration due to gravity [m/s^2]
    mnac = 136000; % mass of the nacelle [kg]
    mrot = 102589; % mass of the rotor [kg]
    mrot = 100850;
    Fzr = mrot.*g; % weight of the rotor [N] 
    Fn = mnac.*g; % force in the z direction due to nacelle weight [N]
    Mfric = Mfric1+Mfric2+Mfric3;    
    M1 = Mxr + Fyr.*zr; % [Nm]
%    M2 = Myr - Fzr*cosd(tilt).*xr + Fxr.*zr - Fn*cosd(tilt).*xn - Fn*sind(tilt).*zn; % working combination
    M2 = Myr - Fzr*cosd(tilt).*xr + (Fxr.*zr.*cosd(tilt) - Fxr.*sind(tilt).*(xr*cosd(tilt))) - Fn.*xn; % trial combination
%     M2 = MyHF.*cosd(tilt).*10^3 + (FzHR.*2.195.*cosd(tilt).*10^3 - FzHR.*(xr*cosd(tilt)).*sind(tilt).*10^3) - (Fn.*xn.*cosd(tilt) + Fn.*zn.*sind(tilt)) + Mfric.*cosd(tilt) - Fzr.*xr.*cosd(tilt); % flex calculation./cosd(tilt)
%    Mtilt = sqrt(M1.^2 + M2.^2);
    figure; % amplitude doesn't match
    plot(s.Timefromstartofoutput, s.YawbearingMy, 'r.-', s.Timefromstartofoutput, M2, 'g.-', Time, MyYB.*10^3, 'b.-');
    legend('Bladed', 'M2', 'Flex');
    title('comparison of yaw bearing tilt moment')
    
    figure;
    plot(s.Timefromstartofoutput, Mtilt, 'k.-', s.Timefromstartofoutput, s.YawbearingMx, 'r.-', s.Timefromstartofoutput, M1, 'g.-')
    legend('Flex', 'Bladed', 'M1');
%    Myaw = Mzr + Fyr*xr + Mbrake + Mfric;
end % Validation

%% Function for bend twist coupling parameters
%------------- END OF CODE --------------
