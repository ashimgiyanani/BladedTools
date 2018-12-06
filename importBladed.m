function [matFile, txtFile, s] = importBladed(varargin)
% Function to import Bladed project run and all the sensors into Matlab

% Input:
%   prjFile: blade project file path
%   Append: Append to txt or mat files
%   AssignVar: Assign the data to variables
%   Format: format of the file to be read. Output of bladed in binary or ascii format
%   version: = 4.6; % option 3.85 or 4.6

% Output:
%   matFile: path to *.mat file location
%   txtFile: path to *.txt file location
  
% Examples:
% clear; clc; close all;
% prjFile = 'c:\Users\ashim.giyanani\OneDrive - Windwise GmbH\Windwise\Bladed\2BEnergy_MC141\MC141_rev47_DoFact_Inflow_min5deg_ecd_p.$PJ'; % bladed project file
% Append = 1; % 1 or 0
% AssignVar = 1; % 1 or 0
% Format = 'binary';    % ascii /binary
% Version = 4.6; % option 3.85 or 4.6


switch nargin
       case 5
            prjFile   = varargin{1};
            Append    = varargin{2};
            AssignVar = varargin{3};
            Format    = varargin{4};
            Version   = varargin{5};
       case 4
            prjFile   = varargin{1};
            Append    = varargin{2};
            AssignVar = varargin{3};
            Format    = varargin{4};
            Version   = 4.6;
       case 3
            prjFile   = varargin{1};
            Append    = varargin{2};
            AssignVar = varargin{3};
            Format    = 'binary'; % ascii /binary
            Version   = 4.6;
       case 2
            prjFile   = varargin{1};
            Append    = varargin{2};
            AssignVar = 1; % 1 or 0
            Format    = 'binary'; % ascii /binary
            Version   = 4.6;
       case 1
            prjFile   = varargin{1};
            Append    = 1;
            AssignVar = 1;
            Format    = 'binary'; % ascii /binary
            Version   = 4.6;
       otherwise
            error('provide a *.prj filepath \n');
end

if isempty(prjFile) || (~isstr(prjFile))
   error('Prj file not found, check the path again \n');
end
pathString = strsplit(prjFile, '\');                                                                     % splitting the file name with the given folder separator '\'
prjString = strsplit(string(pathString(end)), '.');                                                      % getting the project name
folderString = strjoin(string(pathString(1:end-1)),'\');                                                 % combining individual folder elements into a folder path
prjName = prjString(1);                                                                                  % project name
cd(char(folderString));                                                                                        % change current directory
prjInfo = dir(char(strcat(prjName,'.%*')));                                                                    % project info file
prjData = dir(char(strcat(prjName,'.$*')));                                                                    % project data file
idx_remove = [];
for kk = 1:numel(prjInfo)
    infoString(kk,:) = strsplit(prjInfo(kk).name, '%');                                                  % splitting the filename into prj name and sub-module number
    if (strcmp(cellfun(@cellstr,infoString(kk,2)),'AE'))
        idx_remove = [idx_remove,kk];
    elseif (cellfun(@str2num,infoString(kk,2)) == 29) || (cellfun(@str2num,infoString(kk,2)) == 296) ||...
            (cellfun(@str2num,infoString(kk,2)) == 297)
        idx_remove = [idx_remove,kk];
    end
end
prjInfo(idx_remove) = [];
prjData(idx_remove) = [];
clearvars infoString 
CombinedSensors = [];                                                                                    % init
Combined_modNames = {};                                                                                  % init
Combined_varNames = {};                                                                                  % init
Combined_unitNames = {};                                                                                 % init
Combined_vsvNames = {};
for ii = 1:numel(prjInfo)                                                                                % for the number of project info files
    infoString(ii,:) = strsplit(prjInfo(ii).name, '%');                                                  % splitting the filename into prj name and sub-module number
    prjDataFile(ii,:) = strcat(folderString, '\', strjoin(infoString(ii,:), '$'));                       % same sub-module number for data file
    prjInfoFile(ii,:) = strcat(folderString, '\', prjInfo(ii).name);                                     % info file
    % Extract information about the data file from the info file
    fid = fopen(prjInfoFile(ii,:), 'r');                                                                 % open the info file
    for k = 1:8 % number of headers until useful information
      ignoreLine = fgetl(fid);                                                                           % skip lines
    end % no. of headers
    dimLine = fgetl(fid);                                                                                % the number of dimensions given in the info file
    dimSplit = regexp(dimLine,'\s+', 'split');                                                           % split the string
    if size(dimSplit,2) == 3
      NCols = cellfun(@str2double, dimSplit(2));                                                           % no. of columns in the data file
      NRows = cellfun(@str2double, dimSplit(3));                                                           % no. of rows in the data file
      NSubs = 1;
    elseif size(dimSplit,2) == 4
      NCols = cellfun(@str2double, dimSplit(2));                                                           % no. of columns in the data file
      NRows = cellfun(@str2double, dimSplit(3))*cellfun(@str2double, dimSplit(4));                         % no. of rows in the data file, 3 dimensions i.e.
    else
      error('Soemthing wrong there');
    end  
    modLine = fgetl(fid);                                                                                % extract the module line
    modSplit = regexp(modLine,'[\'']+', 'split');                                                        % split a string
    modNames = modSplit(2);                                                                              % name of the sub-module
    varLine = fgetl(fid);                                                                                % get variable list line
    if (Version == 4.6)
        varSplit = regexp(varLine,'[\'']+', 'split');                                                        %
        varNames = varSplit(2:2:end);                                                                        % name of the variables in the data file
    elseif (Version == 3.85)
        varSplit_tab = regexp(varLine,'\s+', 'split');
        varNames = varSplit_tab(2:end);        
        if numel(varNames) > NCols
            has_quotes = regexp(varNames, '[\'']', 'match' );
            idx_quotes = find(~cellfun(@isempty,has_quotes));
            for iq = 1:numel(idx_quotes)/2
                idx_start = idx_quotes(2*iq-1);
                idx_end = idx_quotes(2*iq);
                varNames(idx_start) = cellstr(strjoin(varNames(idx_start:idx_end)));
                varNames(idx_start+1:idx_end) = {''};
            end
        end
    end
    varNames = regexprep(varNames, '[\''\s]+', '');
    varNames(cellfun(@isempty,varNames)) = [];
    
    
    while NCols ~= size(varNames,2)
        varSplit = regexp(varLine, '\s+', 'split');
        varMatch = regexp(varSplit, '[\'']+', 'split');
         idx_join = (cellfun('size',varMatch, 2)~=1);
        varJoin = [varMatch{idx_join}];
        varJoin = varJoin(~cellfun(@isempty,varJoin));
        length_join = FnFind1(idx_join);
        NameSplit = cell(1, numel(length_join));
        for iv = 2:numel(length_join)
            if (length_join(iv) == 0) && (length_join(iv-1) == 0)
               NameSplit(iv) = varSplit(iv);
            elseif length_join(iv) ~= 0
                Nstr = length_join(iv);
                tempvar = strjoin(varSplit(iv:iv+Nstr-1));
                NameSplit(iv) = {regexprep(tempvar, '[\''\s]+', '')};
            elseif (length_join(iv) == 0) && (length_join(iv-1) ~= 0)
                NameSplit(iv) = {[]};
            end
        end
        idx_filled = find(~cellfun(@isempty, NameSplit));
        varNames = NameSplit(idx_filled);
    end
    unitLine = fgetl(fid);                                                                               % get the units line
    unitSplit = regexp(unitLine,'\s+', 'split');                                                         % split a string
    unitNames = unitSplit(2:end);                                                                        % units of the variables in the data file
    unitNames = unitNames(~cellfun(@isempty, unitNames));
    if (contains(prjInfoFile(ii), '.%29')) | (NCols ~= size(unitNames,2))
        disp('Units unreadable, assigning default units');
        unitNames = {'NaN'};
    end
    if size(dimSplit,2) == 4
       ignoreLine = fgetl(fid);
       ignoreLine = fgetl(fid);
       ignoreLine = fgetl(fid);
       subvarLine = fgetl(fid);
       count_multLines = 1;
       subvarNames = [];
       while isempty(regexp(string(subvarLine), 'AXISLAB*', 'match')) || (count_multLines < 2) 
           subvarSplit = regexp(subvarLine,'[\'']+', 'split');
           if size(subvarSplit,2) == 1
                 subvarSplit = regexp(subvarLine,'[\s]+', 'split');
                 subvar = subvarSplit(2:end);
           elseif size(subvarSplit,2) > 1
                 subvar = subvarSplit(2:2:end);   
           end
           subvarNames = [subvarNames, subvar];
           count_multLines = count_multLines + 1;
           subvarLine = fgetl(fid);
           continue
       end
       NSubs = cellfun(@str2double, dimSplit(3));
       vsvName = {};
       for nc = 1:NCols
           for ns = 1:NSubs
              vsvName{ns,nc} = strjoin([varNames(nc), subvarNames(ns)], '_');
           end
       end
       vsvName = reshape(vsvName, [],1);
       varNames = vsvName';
    end
    fclose(fid);                                                                                         % close the file
    clearvars fid

    % Extract the data from the data file
    prjData = cell(NRows, NCols);                                                                        % memory allocation
    formatSpec = repmat('%f', 1, NCols);                                                                 % format of the data file to be extracted
    fid = fopen(prjDataFile(ii,:), 'r');                                                                 % open data file
    if strmatch(Format,'ascii')
      prjData = textscan(fid, formatSpec, 'Delimiter', '\t','EmptyValue', NaN, 'returnOnError', 1);        % read the data file
    elseif strmatch(Format,'binary')
      prjData = fread(fid, NRows*NCols*NSubs, 'float32');
      prjData = num2cell((reshape(prjData, NCols, NRows))');
    end
    fclose(fid);                                                                                         % close the data file
    Combined_modNames = [Combined_modNames, regexprep(modNames, '\s', '')];                              % combine the sub-module names
    Combined_varNames = [Combined_varNames, regexprep(varNames, '\s', '')];                              % combine the variable names
    Combined_unitNames = [Combined_unitNames, unitNames];                                                % combine the units for the variables
%   Combined_vsvNames = [Combined_vsvNames, regexprep(vsvName, '\s', '')];
    if size(dimSplit,2) == 4
        new_prjData = [];
        for ic = 1:NCols
            temp_prjData = (reshape(cell2mat(prjData(:,ic)), NSubs, NRows/NSubs))';
            new_prjData = [new_prjData, temp_prjData];
        end
        prjData = new_prjData;
        clearvars new_prjData temp_prjData
    elseif size(dimSplit,2) == 3
        prjData = reshape(cell2mat(prjData),NRows/NSubs, NCols*NSubs);
    end
    CombinedSensors = [CombinedSensors, prjData];                                                     % Combine the data from the sensors
    clearvars prjData
    fprintf('[%s]:Extracted %d/%d files \n', datetime, ii, numel(prjInfo));
end % no. of prjInfo files

if Append == 1
    txtFile = sprintf('combined_%s.txt',prjName);
    fid = fopen(txtFile,'a');
    headerString = [strjoin(Combined_varNames, ';') '\n'];
    formatSpec = [repmat('%.2f;',1,size(CombinedSensors,2)-1), '%.2f\r\n'];
    fprintf(fid, headerString);
    fprintf(fid, formatSpec, CombinedSensors);
    fclose(fid);
    fprintf('[%s]:Saved to *.txt file (%s)\n', datetime, txtFile);
end

if AssignVar == 1
    matFile = sprintf('combined_%s.mat',prjName);
    save(matFile);
    Combined_varNames = regexprep(Combined_varNames, '\.', 'p');
    Combined_varNames = regexprep(Combined_varNames', '-|\(', '_');
    Combined_varNames = regexprep(Combined_varNames', '\)', ''); 
    for ival = 1:numel(Combined_varNames)
          s.(char(Combined_varNames(ival))) = CombinedSensors(:,ival);
    end
    save(matFile);
    fprintf('[%s]:Saved to *.mat file (%s)\n', datetime, matFile);
end

%%%***************Extended functions********************%%%%%%%%%%%%%%%
function out = FnFind1(a)
  out     = zeros(size(a));
  aa      = [0,a,0];
  idx      = strfind(aa, [0 1]);
  out(idx) = strfind(aa, [1 0]) - idx;
end% function

% function FnCreateVars(Combined_varNames,CombinedData)
% % script to create variables and exchange between two workspaces
% evalin('base', [char(Combined_varNames(ival)) '= CombinedSensors(:,ival);']);
% end

end % function
