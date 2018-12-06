function [cum_locTwist, bladed_meanSig, bladed_Psi, bladed_Tpsi, flex_meanSig, flex_Psi, flex_Tpsi, flex_r, bladed_r] = FnBendTwist(ExcelFilename, matFile, Time, Psi, queryPsi, queryRotorazimuthangle, stepStart, currentMx, currentFz)
%FnImportFile - to calculate the bend-twist coupling based on the profile details
%
%
% Syntax:  [cum_locTwist, bladed_meanSig, flex_meanSig] = FnBendTwist(ExcelFilename, matFile, Time, Psi, queryPsi, queryRotorazimuthangle, currentMx, currentFz)
%
% Inputs:
%    ExcelFilename - e.g. 'p:\02_Projekte\08_Maxcap\04_Rotorblade\Sinoi_LZ69-3.0\LZ69-3.0-Decrypt-Excl.Ice-Rev10.11-GH4.3-20170726_CBe.xlsx' retain format as on 27/11/2018
%    matFile - link to *.mat file 'c:\Users\ashim.giyanani\OneDrive - Windwise GmbH\Windwise\Bladed\2BEnergy_MC141\combined_MC141_rev46_DoFact_Inflow_min5deg_ecd_n.mat'
%    Time - input time from Flex [s]
%    Psi  - rotor position in deg from Flex simulation [deg]
%    queryPsi  - query rotor position to match the rotor position of Bladed i.e. done to avoid variations due to different azimuth positions, Flex Blade1 faces up at 180 deg [deg]
%    queryRotorazimuthangle - query rotor position so that the aimuth angle for Flex and Bladed is the same [rad]
%    currentMx - just to exchange information between the matlab workspace and function's local workspace, initialize as currentMx=[]
%    currentFz - same as currentMx, initialize as currentFz=[]
%    stepStart - begnning time of the steo in the time series
%
% Outputs:
%    cum_locTwist - Twist time series for all radial locations in format [radial pos, time] [deg]
%    bladed_meanSig - mean Twist at a radial position for the query Rotorazimuth angle [deg]
%    flex_meanSig - mean Twist at a radial position for the query Psi (rotor azimuth angle) [deg]
%    flex_Psi - Rotor azimuth position matching the criteria and is close to bladed_Psi [deg]
%    bladed_Psi - Rotor azimuth position matching the criteria and is close to flex_Psi [deg]
%    flex_r -  radial positions from flex [m]
%    bladed_r -  radial positions from bladed [m]
%    
% Example:
%    Line 1 of example
%    Line 2 of example
%    Line 3 of example
%
% Other m-files required: CompareSignals.m, CompareSignals_Flex.m, CompareSignals_Bladed.m, CompareSignals_BendTwist.m 
% Subfunctions: none
% MAT-files required: matFile from Bladed output, ExcelFileName for blade configuration, 
%
% See also: CompareSignals.m, CompareSignals_Flex.m, CompareSignals_Bladed.m, CompareSignals_BendTwist.m 

% See also: other function name, other functions
% Author: Ashim Giyanani, WindWise gmbh, Muenster
% WindWise gmbh, Muenster
% Created by      : Ashim Giyanani, PhD, Lidar applications Wind Energy, TU Delft, the Netherlands
% Email           : ashimgiyanani@yahoo.com
% Github          : 
% Created         : 27.11.2018 18:15:32
% Last Updated    : 27.11.2018 18:15:39
% Updated by      : Ashim Giyanani
% Comments        : Final v1.0

%------------- BEGIN CODE --------------


% Script 
  % import the airfoil details
%   filename = 'p:\02_Projekte\08_Maxcap\04_Rotorblade\Sinoi_LZ69-3.0\LZ69-3.0-Decrypt-Excl.Ice-Rev10.11-GH4.3-20170726_CBe.xlsx';
  filename = ExcelFilename;  
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
  flex_r = ndata(:,1);    
  
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
  
  % find the azimuth position where the rotor position is at 0 deg or where the rotor position matches criteria for Bladed rotor as well
  condn1 = (round(Psi) == round(queryPsi));               % find rotor position meeting the query angle
  condn2 = (Time > (stepStart-2) & Time < (stepStart+2)); % find the position within a specified range
  [temp_idx1,~] = find(condn1 & condn2);                  % both conditions
  [~, idx] = min(abs(stepStart - Time(temp_idx1)));       % find the closest position to the step begin time
  idx1_azimuthpos = temp_idx1(idx);                       % assignment
  nPos = 41;
  cum_Mx1=[];
  cum_Fz1 = [];
  for in = 1:nPos
%      evalin('caller', ['currentMx =' sprintf('Mx1%d', in);]);
      currentMx = evalin('base', sprintf('Mx1%d', in));
      cum_Mx1 = [cum_Mx1, currentMx];
%      evalin('caller', ['currentFz =' sprintf('Fz1%d', in);]);
      currentFz = evalin('base', sprintf('Fz1%d', in));
      cum_Fz1 = [cum_Fz1, currentFz];
  end
  clearvars condn1 condn2 temp_idx1 idx 
  
  Mx = cum_Mx1' + cum_Fz1'.*y_shear_centre; % moment around the shear centre for a radial location or for a particular time
  dr = ([diff(Length);0])./2+([0;diff(Length)])./2;  % blade element length
  local_Twist = rad2deg(dr.*Mx.*10^3./Tor_S); % local twist in the element
  cum_locTwist = cumsum(local_Twist); % cumulative twist in the blade until the tip
  
  % Compute the radial twist along the blade
    % bladed
  String = 'Blade1rotationaboutz_plane'; % Blade1x_deflection_perpendiculartorotorplane, Blade1rotationaboutz_plane
  s1 = load(matFile, 's');
  fnames1 = fieldnames(s1.s);
  filt_s1 = rmfield(s1.s, fnames1(find(cellfun(@isempty,strfind(fnames1, String)))));
  fnames1 = fieldnames(filt_s1);
  temp_str1 = regexp(fnames1, '\d+p\d+', 'match');
  for i = 1:size(temp_str1,1)
      temp_str1{i} = regexprep(temp_str1{i}, 'p', '.');
      bladed_r(i) = cellfun(@(x) (str2num(x)), temp_str1{i}, 'UniformOutput', true);
  end
  bladed_r = double(bladed_r);
  
  unit = 'rad'; % option of 'rad' or 'deg' or 'm'
  for ii = 1:numel(fnames1)
        if strmatch(unit, 'rad')
           bladed_signal = rad2deg(filt_s1.(fnames1{ii}));
        end
        condn1 =  (round(rad2deg(s1.s.Rotorazimuthangle)) == round(queryRotorazimuthangle));
        condn2 =  (s1.s.Timefromstartofoutput > (stepStart-2) & s1.s.Timefromstartofoutput < (stepStart+2) );
        [temp_idx2, ~] = find((condn1 & condn2));
        [~, idx] = min(abs(stepStart - s1.s.Timefromstartofoutput(temp_idx2)));
        idx2_azimuthpos = temp_idx2(idx);
        bladed_meanSig(ii) = nanmean(bladed_signal(idx2_azimuthpos));
        temppos{ii} = regexp(fnames1{ii}, '\d*', 'match');
        radpos(ii) = cellfun(@str2double, temppos{ii}(2));
  end
     % Flex
  flex_meanSig = nanmean(cum_locTwist(:,idx1_azimuthpos),2);
  flex_Psi = Psi(idx1_azimuthpos);
  flex_Tpsi = Time(idx1_azimuthpos);
  bladed_Psi = s1.s.Rotorazimuthangle(idx2_azimuthpos);
  bladed_Tpsi = s1.s.Timefromstartofoutput(idx2_azimuthpos);

  flex_signal = cum_locTwist(6,:); %#####
  bladed_time = s1.s.Timefromstartofoutput;
  flex_time = Time;

