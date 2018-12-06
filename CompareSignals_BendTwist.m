% Script to compare the bend twist coupling results from Flex and bladed
% Other m-files required: CompareSignals.m, compareSignals_Bladed, CompareSignals_Flex
% Subfunctions: none
% MAT-files required: none
%% Files for Bend twist coupling tests
% bladedFile1 = 'w:\Bladed\2Benergy_MC141\combined_MC141_rev44_WG_061118_DoFactive_NoVar6p3.mat';
% bladedFile2 = 'w:\Bladed\2Benergy_MC141\combined_MC141_rev44_WG_061118_DoFdeactivated_NoVar6p3.mat';
% FlexFile1 = 'c:\Users\ashim.giyanani\OneDrive - Windwise GmbH\Windwise\Flex\ExtremeWindShear_rev04_126mHH\Positive_EWM_Inflow_min5deg.txt';
% FlexFile2 = 'c:\Users\ashim.giyanani\OneDrive - Windwise GmbH\Windwise\Flex\ExtremeWindShear_rev04_126mHH\Positive_EWM_DOF_deactivated_Inflow_min5deg.txt';
%
% See also: OTHER_FUNCTION_NAME1,  OTHER_FUNCTION_NAME2

% Author: Ashim Giyanani
% Windwise GmbH, Münster
% email address: ashimgiyanani@gmail.com
% Website: n/a
% Nov 2018; Last revision: 20.11.2018 14:27:52

% Steps to be followed
% #a1 Import Bladed and Flex output to be compared
% #b1 plot the plots of pitch angle, wind speed, electrical power, rotor speed, thrust force, thrust force - not tilted, tower bottom bending moment, rotor torque,

%------------- BEGIN CODE --------------

%% Files for Bend twist coupling tests
% bladedFile1 = 'w:\Bladed\2Benergy_MC141\combined_MC141_rev44_WG_061118_DoFactive_NoVar6p3.mat';
% bladedFile2 = 'w:\Bladed\2Benergy_MC141\combined_MC141_rev44_WG_061118_DoFdeactivated_NoVar6p3.mat';
% FlexFile1 = 'c:\Users\ashim.giyanani\OneDrive - Windwise GmbH\Windwise\Flex\ExtremeWindShear_rev04_126mHH\Positive_EWM_Inflow_min5deg.txt';
% FlexFile2 = 'c:\Users\ashim.giyanani\OneDrive - Windwise GmbH\Windwise\Flex\ExtremeWindShear_rev04_126mHH\Positive_EWM_DOF_deactivated_Inflow_min5deg.txt';

clear; clc; close all;
addpath('c:\Users\ashim.giyanani\OneDrive - Windwise GmbH\Windwise\Matlab');
PlotFig = 0;
SaveFig = 0;
Validation = 0; 
impBladed = 0; % import bladed into mat file and the workspace 0 - no import, 1 - full import
oldFile = 0; % old imports of the Bladed format were not saved to the struct variable  's'
SpecFig = 1;

% Prepare the Bladed files
if impBladed == 1
  Append = 1;
  AssignVar = 1;
  Format= 'binary';
  Version = 4.6;
  prjFile = 'w:\Bladed\2Benergy_MC141\MC141_rev44_WG_061118_DoFactive_NoVar6p3.$PJ';
  [matFile, ~, s] = importBladed(prjFile, Append, AssignVar, Format, Version);
else
  matFile = 'w:\Bladed\2Benergy_MC141\combined_MC141_rev44_WG_061118_DoFactive_NoVar6p3.mat';
  load(matFile)
end

% Input
Delimiter = '\t'; % delimiter in the file to be read
nHead = 14; % number of lines before header names
fpath = 'c:\Users\ashim.giyanani\OneDrive - Windwise GmbH\Windwise\Flex\ExtremeWindShear_rev04_126mHH\Positive_EWM_Inflow_min5deg.txt';
String = 'VuHc|Teta2|OmRot|P_Gen|FzHR|FzYB|MyTB|MzHR|FzH2|MyYB|MyHR|MxHF|TwrClr|Mx*|My*|Mz*|Fx*|Mfric*|\wH1$|M*\w11$';
String = '*Tip*|Uf*|Uk*|TwrClr|FTip*|Mx*|Fz*|Psi';

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

% import the airfoil details
filename = 'p:\02_Projekte\08_Maxcap\04_Rotorblade\Sinoi_LZ69-3.0\LZ69-3.0-Decrypt-Excl.Ice-Rev10.11-GH4.3-20170726_CBe.xlsx';
fid  = fopen(filename, 'r');
[ndata, textdata, alldata] = xlsread(filename, 6);
fclose(fid);
ndata = [ndata(1,:);ndata];
% replace certain column names that do not pass the criteria
invalidColNames = {'t25 centre relative to local axis [%]',...
                        'ratio of distance from t25 centre to pitch axis and chord length [-]',...
                               'y_shear_centre  [m]', ...
                                   'y_t25_centre - y_shear_centre  [m]',...
                                                 'd(y_t25_centre, y_shear_centre)  [rel]'};
replaceColNames = {'AC_t25_local', 'Ratio_t25_p_c', 'CS_yp', 'dy_AC_CSyp', 'Ratio_dy_AC_CSyp_chord'};
for ic = 1:numel(invalidColNames)
    idx_invalid = strcmp(textdata,invalidColNames{ic});
    textdata{idx_invalid} = replaceColNames{ic};
end    

% Arrange the names in a way such that variable names in Matlab can be created
tempNames = regexp(textdata, '[|/|\s*', 'split');
for in = 1:size(tempNames,2)
    if in == 8
       continue;
    end
    fprintf('Variables names assignment')
    fprintf('%s -> %s \n',string(textdata{in}),string(tempNames{in}(1)));
    eval([char(tempNames{in}(1)) '= ndata(:,in);']);
end
fprintf('[%s] - Airfoil import completed \n', datetime);

CE_Rxy = sqrt(CE_X.^2 + CE_Y.^2); % shortest distance between the local axis [%]
relAC = CE_Rxy - 25; % relative distance between the aerodynamic centre (at c/4 location) [%]
Ratio_ACpitch_chord = max((-Ny./Chord.*100+ relAC),-100)./100; % ratio of distance from t25 centre to pitch axis and chord length [-]
y_shear_centre = Ny + Chord.*(CS_Y - CE_Y)./100; % distance of shear centre from neutral axis
dy_CS_AC = (-relAC.*Chord) - y_shear_centre; % relative distance between the Aerodynamic centre and shear centre
rel_dy_CS_AC = dy_CS_AC./Chord; % relative distance to chord ratio

% find the azimuth position where the rotor position is at 0 deg
[idx1_azimuthpos_0deg, val1_azimuthpos_0deg] =find(round(Psi,-1) == 0);
idx_Time = find(Time == 60.16)
nPos = 41;
cum_Mx1=[];
cum_Fz1 = [];
for in = 1:nPos
    eval(['currentMx =' sprintf('Mx1%d;', in);]);
    cum_Mx1 = [cum_Mx1, currentMx];
    eval(['currentFz =' sprintf('Fz1%d;', in);]);
    cum_Fz1 = [cum_Fz1, currentFz];
end

% CAlculate the twist is ever blade section from Flex signals
Mx = cum_Mx1' + cum_Fz1'.*y_shear_centre; % moment around the shear centre for a radial location or for a particular time
dr = ([diff(Length);0])./2+([0;diff(Length)])./2;  % blade element length
local_Twist = rad2deg(dr.*Mx.*10^3./Tor_S); % local twist in the element
cum_locTwist = cumsum(local_Twist); % cumulative twist in the blade until the tip

%% Calculate the twist in each blade section from Bladed signals
% sort out the required signal from the struct
StringM = 'Blade1Mz_Rootaxes';
StringF = 'Blade1Fx_Rootaxes';
s2 = load(matFile,'s');
fnames2 = fieldnames(s2.s);
filtM_s2 = rmfield(s2.s, fnames2(find(cellfun(@isempty,strfind(fnames2, StringM)))));
filtF_s2 = rmfield(s2.s, fnames2(find(cellfun(@isempty,strfind(fnames2, StringF)))));
fnamesM = fieldnames(filtM_s2);
fnamesF = fieldnames(filtF_s2);
Bld_Mz1 = [];
Bld_Fx1 = [];
unit = 'rad';
for ii = 1:numel(fnamesM)
       Bld_Mz1 = [Bld_Mz1, filtM_s2.(fnamesM{ii})];
       Bld_Fx1 = [Bld_Fx1, filtF_s2.(fnamesF{ii})];
end
[idx2_azimuthpos_0deg, val2_azimuthpos_0deg] = find(round(rad2deg(s.Rotorazimuthangle),-1) == 0);
Bld_Mz = -Bld_Mz1' + Bld_Fx1'.*y_shear_centre; % moment around shear centre calculated from excel sheet, alternative to use the ones defined in Bladed
Bld_local_Twist = rad2deg(dr.*Bld_Mz./Tor_S); % local twist in the element
Bld_locTwist = cumsum(Bld_local_Twist); % cumulative twist along the blade length in Bladed 
Bld_meanSig = nanmean(Bld_locTwist(:,idx2_azimuthpos_0deg),2);

% Compute the radial twist along the blade
% bladed
String = 'Blade1rotationaboutz_plane'; % Blade1x_deflection_perpendiculartorotorplane, Blade1rotationaboutz_plane
s1 = load(matFile, 's');
fnames1 = fieldnames(s1.s);
filt_s1 = rmfield(s1.s, fnames1(find(cellfun(@isempty,strfind(fnames1, String)))));
fnames1 = fieldnames(filt_s1);
unit = 'rad'; % option of 'rad' or 'deg' or 'm'
for ii = 1:numel(fnames1)
      if strmatch(unit, 'rad')
         bladed_signal = rad2deg(filt_s1.(fnames1{ii}));
      end
      [idx2_azimuthpos_0deg, val2_azimuthpos_0deg] = find(round(rad2deg(s1.s.Rotorazimuthangle),-1) == 0);
      bladed_meanSig(ii) = nanmean(bladed_signal(idx2_azimuthpos_0deg));
      temppos{ii} = regexp(fnames1{ii}, '\d*', 'match');
      radpos(ii) = cellfun(@str2double, temppos{ii}(2));
end
% Flex
flex_meanSig = nanmean(cum_locTwist(:,idx1_azimuthpos_0deg),2);

flex_signal = cum_locTwist(6,:); %#####
bladed_time = s.Timefromstartofoutput;
flex_time = Time;

% plotting the time series of the signal 1 (bladed) vs signal 2 (Flex)
p1 = figure;
plot(Length, Bld_meanSig, 'm-*', radpos, bladed_meanSig, 'r.-', Length, flex_meanSig, 'k.-');
title('Twist about z axis for rotor position = 180 deg'); % #####
xlabel('Radial position [m]');
ylabel('Twist about z-axis [deg]');
legend('Bld_calc', 'Bladed', 'Flex');
grid minor
