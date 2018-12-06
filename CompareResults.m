% Script to scan through the folder and pick up the signal of importance

clear; clc; close all
addpath('c:\Users\ashim.giyanani\OneDrive - Windwise GmbH\Windwise\Matlab');
% list all the files required
PlotFig = 1;
SaveFig = 1;

% Get the Flex results
Delimiter = '\t'; % delimiter in the file to be read
nHead = 14; % number of lines before header names
fpath = 'c:\Users\ashim.giyanani\OneDrive - Windwise GmbH\Windwise\Flex\StairCase_rev02\Treppe_AllDOFdeactivated_GenDelay0p1s_cd0_rev03.txt';
String = 'VuHc|Teta2|OmRot|P_Gen|FzHR|FzYB|MyTB|MzHR|FzH2|MyYB|MyHR|MxHF|Mx*|My*|Mz*|\wH1$|M*\w11$';
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

% get the Bladed results
fpath = 'w:\Bladed\2Benergy_MC141\';
files = dir(fullfile(fpath,'*.mat'));
MatFiles = strcat(fpath,{files.name}');

for i = 1:numel(files)
  filename = MatFiles{i};
  tempSplit = regexp(filename, '_', 'split');
  prjName = regexp(tempSplit{end}, '.mat', 'split');
  prjName = string(prjName(1));
  load(filename)%,...
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
                      MzTB];                %25 foundations roll tower base
  
      PlotVar_Bladed = [
                        Hublongitudinalwindspeed,...                       %1 wind speed                                  -> matches          
                        Rotorspeed.*30./pi,...                             %2 rotor speed                                 -> matches        
                        Electricalpower./10^3,...                          %3 electrical power generated                  -> matches  
                        rad2deg(Blade2pitchangle),...                      %4 pitch angle                                 -> matches       
                        StationaryhubFx./10^3,...                          %5 thrust                                      -> matches      
                        StationaryhubMx./10^3,...                          %6 torque                                      -> matches 
                        FoundationMy./10^3,...                             %7 tower base tilt moment                      -> matches (NCoMmin1p6m)
                        YawbearingFx./10^3,...                             %8  yaw thrust                                 -> matches
                        Bladeroot1Fx./10^3,...                             %9 blade oop bending force                     -> matches
                        YawbearingMy./10^3,...                             %10 yaw bearing tilt moment                    -> ### differs a bit, offset 50-100 kNm at below-rated 
                        RotatinghubMy./10^3,...                            %11 rotating hub bending moment                -> matches
                        -StationaryhubMz./10^3,...                         %12 fixed hub yaw moment                       -> ### differs a bit offset 100kNm
                        Bladeroot1Mx./10^3,...                     %13 in plane bending moment                    -> matches 
                        -Bladeroot1My./10^3,...                    %14 oop bending moment                         -> matches 
                        Bladeroot1Mz./10^3,...                     %15  torsional moment                          -> ### differs significantly
                        Bladeroot1Fx./10^3,...                             %16 in plane bending force                     -> matches  
                        -Bladeroot1Fy./10^3,...                            %17 out of plane bending force                 -> matches    
                        Bladeroot1Fz./10^3,...                             %18 torsional force                            -> matches   
                        RotatinghubMz./10^3,...                            %19 bending moment                             -> matches        
                        RotatinghubMx./10^3,...                            %20 torque                                     -> matches               
                        StationaryhubMy./10^3,...                          %21 fixed hub tilt moment                       -> ### big offset around 100 kNm
                        -YawbearingMz./10^3,...                            %22 yaw bearing yaw moment                     -> ### differs a bit     
                        YawbearingMx./10^3,...                             %23 roll                                       -> matches                   
                        -FoundationMz./1000,...                            %24 yaw tower base                             -> matches             
                        FoundationMx./1000,...                             %25 roll tower base                            -> ### differs a bit             
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
                   'FoundationMx'};                                        %25 foundations roll tower base
  
      plotYlabels = {...
                     'hub height wind speed (m/s)',...
                     'rotor speed (rpm)',...
                     'generator power (kW)',...
                     'pitch angle (deg)',...
                     'thrust force (kN)',...
                     'Torque (kNm)',...
                     'Tower base tilt moment (kNm)',...
                     'Thrust, yaw bearing (kN)', ...
                     'blade oop fore (kN)',...
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
                     'FoundationMx'};
  
  
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
                   'FoundationMx'};
                   
      dt = mode(diff(Time));
      Fs = 1/dt;
      selVars = [10,15,21]; % or 1:25 referring to the variables
      for ip = 1:size(selVars,2)
          sv = selVars(ip);
          p1 = figure();
          plot(Timefromstartofoutput, PlotVar_Bladed(:,sv), 'r.-',Time,  PlotVar_Flex(:,sv), 'k.-');
          title(plotTitle{sv});
          xlabel('Time (s)');
          ylabel(plotYlabels{sv});
          legend('Bladed','Flex');
          grid minor
          if SaveFig == 1
             saveFolder = 'c:\Users\ashim.giyanani\OneDrive - Windwise GmbH\Windwise\Reports\Pics\';
             saveFile = sprintf('%sRes_%s_%s',saveFolder, prjName, fileNames{sv});
             saveas(p1, saveFile,  'epsc');
  %          export_fig(p1, saveFile,'-transparent','-q10', '-depsc');
  %            cleanfigure;
  %            matlab2tikz('interpretTickLabelsAsTex',true,'maxChunkLength',15000,'floatFormat', '%.4g', 'checkForUpdates', true,'imagesAsPng',true,  [saveFile,'.tikz'])
          end % savefig       
          
  %         p2 = figure();
  %         [Pf, ff] = periodogram(PlotVar_Flex(:,sv), [], length(PlotVar_Flex(:,sv)), Fs, 'onesided', 'psd');
  %         [Pb, fb] = periodogram(PlotVar_Bladed(:,sv), [], length(PlotVar_Bladed(:,sv)), Fs, 'onesided', 'psd');
  %         loglog(ff, Pf, 'k.-', fb, Pb,'r.-');
  %         legend('Flex', 'Bladed');
  %         ylabel(sprintf('%s Spectra', plotYlabels{sv}),'FontSize', 12);
  %         xlabel('Frequency Hz', 'FontSize', 12);
  %         title(sprintf('%s Spectrum comparison', plotTitle{sv}));
  %         grid on; 
  %         if SaveFig == 1
  %           saveas(p2, sprintf('spect_%s_%s',prjName, fileNames{sv}), 'epsc')
  %         end % savefig
      end % number of plots
  end % plotfig
  fprintf('[%s]: file processes %s \n', datetime, prjName);
  prompt = 'Do you want to move to next file? Y/N [Y]: \n';
  str = input(prompt,'s');
  if strcmpi(str,'y')
     proceed = 1;
     continue;
  elseif strcmpi(str,'n')
         proceed = 0;
         error('Manual stop by pressing no');
  end
  close all
end % number of files

%------------- END OF CODE --------------


