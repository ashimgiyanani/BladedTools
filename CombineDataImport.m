%function [output1,output2] = function_name(input1,input2,input3)
%FUNCTION_NAME - Imports data from the whole dataset and writes the required variables in a *.mat file
%Optional file header info (to give more details about the function than in the H1 line)
%Optional file header info (to give more details about the function than in the H1 line)
%Optional file header info (to give more details about the function than in the H1 line)
%
% Syntax:  [output1,output2] = function_name(input1,input2,input3)
%
% Inputs:
%    input1 - Description
%    input2 - Description
%    input3 - Description
%
% Outputs:
%    output1 - Description
%    output2 - Description
%
% Example:
%    Line 1 of example
%    Line 2 of example
%    Line 3 of example
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: OTHER_FUNCTION_NAME1,  OTHER_FUNCTION_NAME2

% Author: Ashim Giyanani, Ph.D. candidate, Wind Energy
% Wind energy department, TU Delft
% email address: ashimgiyanani@gmail.com
% Website: n/a
% April 2018; Last revision: 20-04-2018

clear; close all; clc;
addpath('d:\Main\MATLAB\Functions\')

tic
%%Providing Default data
Delimiter = ';';
StartRow = 4;
DefinitePath = 0;
DeffPath = 'l:\awep\we\Ashim\Main\Data\TUDelft_XD115_signals_64Hz\2013\month 11\data4Hz.2013-11-31 23.50.00.csv';
TableData = 0;
Append = 1;
ResampleData = 1;
dt64 = 1/64;
dt4 = 1/4;
T = 600;
t64 = (dt64:dt64:T)';
t4 = (dt4:dt4:T)';
count_Append = 0;
InterpolateNaN = 1;
files = {};
nanfiles1 = {};
nanfiles2 = {};

% Get the filenames for XD115 turbine signals
fpath = 'l:\awep\we\Ashim\Main\Data\TUDelft_XD115_signals_64Hz';
if DefinitePath == 0
  % adding the paths for the files to the matrix with name 'files'
  [sub, fls] = subdir(fpath);
  % files = dir(fullfile(fpath,'data4Hz.201*.csv'));
  for i = 1:length(sub)
        for j = 1: length(fls{i})
          files(j,i) = (strcat(sub(i),'\',fls{i}{j}));                                             % all the files in the directory
        end
  end
  files = reshape(files, [],1);
  files1 = files(~cellfun('isempty', files));
else
  files1 = cellstr(DeffPath);
end
clearvars sub fls files i j

% Providing default data for files additional signals file from XD115
fpath = 'l:\awep\we\Ashim\Main\Data\TUDelft_XD115_additional signals_64Hz';
if DefinitePath == 0
  % adding the paths for the files to the matrix with name 'files'
  [sub2, fls2] = subdir(fpath);
  % files = dir(fullfile(fpath,'data4Hz.201*.csv'));
  for i = 1:length(sub2)
        for j = 1: length(fls2{i})
          files(j,i) = (strcat(sub2(i),'\',fls2{i}{j}));                                             % all the files in the directory
        end
  end
  files = reshape(files, [],1);
  files2 = files(~cellfun('isempty', files));
else
  files2 = cellstr(DeffPath);
end
clearvars sub fls files i j

% Providing default data for files Lidar Avent Lidar
fpath = 'd:\Main\Data\data4Hz_WI_AvPr';
String3 = '^Wt10\w*D\w*(Speed$|RWS$|RWS0$)';
if DefinitePath == 0
  % adding the paths for the files to the matrix with name 'files'
  [sub3, fls3] = subdir(fpath);
  % files = dir(fullfile(fpath,'data4Hz.201*.csv'));
  for i = 1:length(sub3)
        for j = 1: length(fls3{i})
          files(j,i) = (strcat(sub3(i),'\',fls3{i}{j}));                                             % all the files in the directory
        end
  end
  files = reshape(files, [],1);
  files3 = files(~cellfun('isempty', files));
else
  files3 = cellstr(DeffPath);
end
clearvars sub fls files i j

%% Creating a default Variable to inlcude the column headers as tested manually
NCols = 13;
FileId = fopen(char(files1(1)), 'r');
DefVarNames1 = textscan(FileId, '%s', NCols, 'delimiter', ';', 'headerlines', 1);
fclose(FileId);
% find the indexes of the colum to be imported
String = '^Wt10\w*(pitch|NaWs|Pel|rotorspeed|RootMoment|AeroTorque|PgspdSse|Vw)\w*';
Var = [DefVarNames1{:}];
ind = regexp(Var, String, 'match');                                                              % regexp to match the variable names to the string input
indf = find(not(cellfun('isempty', ind)));                                                       % find matching instances
index = [indf];                                                                                % indexes to be importedd
VarNames1 = Var(index);                                                                          % variables that will be imported
indfs = indf-1;                                                                                  %
fmt1 = {'%*s'};                                                                                   %
fmt2 = repmat({'%*f'},1, NCols-1);                                                               %
fmt2(indfs) = {'%f'};                                                                            %
FormatSpec1 = [fmt1{:} fmt2{:} '%[^\r\n]'];                                                      % format specification for the import variables
clearvars ind indf index indfs fmt1 fmt2 FileId NCols sub fls files i j

% find the indexes of the colum to be imported
NCols = 12;
FileId = fopen(char(files2(1)), 'r');
DefVarNames2 = textscan(FileId, '%s', NCols, 'delimiter', ';', 'headerlines', 1);
fclose(FileId);
String = '^Wt10\w*(pitch|NaWs|Pel|rotorspeed|RootMoment|AeroTorque|PgspdSse|Vw)\w*';
Var = [DefVarNames2{:}];
ind = regexp(Var, String, 'match');                                                              % regexp to match the variable names to the string input
indf = find(not(cellfun('isempty', ind)));                                                       % find matching instances
index = [indf];                                                                                % indexes to be importedd
VarNames2 = Var(index);                                                                          % variables that will be imported
indfs = indf-1;                                                                                  %
fmt1 = {'%*s'};                                                                                   %
fmt2 = repmat({'%*f'},1, NCols-1);                                                               %
fmt2(indfs) = {'%f'};                                                                            %
FormatSpec2 = [fmt1{:} fmt2{:} '%[^\r\n]'];                                                      % format specification for the import variables
clearvars ind indf index indfs fmt1 fmt2 FileId NCols sub2 fls2 files i j

% find the indexes of the colum to be imported
NCols = 287;
FileId = fopen(char(files3(1)), 'r');
DefVarNames3 = textscan(FileId, '%s', NCols, 'delimiter', ';', 'headerlines', 1);
fclose(FileId);
Var = [DefVarNames3{:}];
ind = regexp(Var, String3, 'match');                                                              % regexp to match the variable names to the string input
indf = find(not(cellfun('isempty', ind)));                                                       % find matching instances
index = [1;indf];                                                                                % indexes to be importedd
VarNames3 = Var(index);                                                                          % variables that will be imported
indfs = indf-1;                                                                                  %
fmt1 = {'%s'};                                                                                   %
fmt2 = repmat({'%*f'},1, NCols-1);                                                               %
fmt2(indfs) = {'%f'};                                                                            %
FormatSpec3 = [fmt1{:} fmt2{:} '%[^\r\n]'];                                                      % format specification for the import variables
clearvars ind indf index indfs fmt1 fmt2 FileId NCols sub2 fls2 files i j


%% Data import from each file in the directory
count = 0;
Teff_series = [];
NewData = [];
Data = [];
nfiles = numel(files1);                                                                          % no. of files to be imported
% nfiles = 2;                                                                                    % manual entry
                                                                      % no. of variables to be imported
for ii = 1956:nfiles % no. of files to be imported
    LoopData = nan(2400,59);
    [fid, ~] = fopen(files1{ii},'r');
    NewData = textscan(fid, FormatSpec1, 'Delimiter', ';', 'EmptyValue', NaN,...
               'headerlines', StartRow-1, 'returnonError',0, 'EndOfLine', '\r\n');                 % import data
    fclose(fid);
    cellsize = cellfun(@numel, NewData(1:end-1));
    if ~all(cellsize==mode(cellsize));
        idx = find(cellsize(:)~=mode(cellsize));
        for jc = 1:numel(idx)
            temp = nan(mode(cellsize),1);
            temp(1:cellsize(idx(jc))) = NewData{idx(jc)};
            NewData{idx(jc)} = temp;
            clearvars temp
        end
    end
    % Extract the date string for verification
    folderstr1 = strsplit(string(files1(ii)),'\');
    filestr1 = strsplit(folderstr1(end), 'Hz.');
    datestring1 = filestr1(end);
%     ts = datenum([NewData{:,1}]);                                                                  % time stamp for wind speed
    Data = cell2mat(NewData(1:end-1));                                                        % colelcting in one array for writing purposes
%    clearvars NewData
    % Check for NaNs and resample the columns
    if InterpolateNaN == 1
        if (sum(sum(isnan(Data))) > 0) && isempty(find(sum(isnan(Data)) > 0.9*numel(t64)))
        col_nan = find(sum(isnan(Data)) > 0);
            for ic = 1:numel(col_nan)
                method='linear';
                idx_nan = isnan(Data(:,col_nan(ic)));
                if idx_nan(1)==1
                    Data(1,col_nan(ic)) = nanmean(Data(1:128,col_nan(ic)));
                elseif idx_nan(end)==1
                    Data(end,col_nan(ic)) = nanmean(Data(end-128:end,col_nan(ic)));
                end
                idx_nan = isnan(Data(:,col_nan(ic)));
                Data(idx_nan,col_nan(ic)) = interp1(t64(~idx_nan), Data(~idx_nan,col_nan(ic)), t64(idx_nan), method, 'extrap');
            end % ic #replace nans condn
        elseif ~isempty(find(sum(isnan(Data)) > 0.9*numel(t64)))
            disp('Too many Nans');
            nanfiles1 = [nanfiles1; files1(ii)];
        end
    end
    
    % Assign the data columns to the variable name after successful interpolation and resampling
    nVar1 = length(VarNames1);
    totalVar = nVar1;
    [nr,nc] = size(Data);
    ResampledData = nan(round(nr/16)+8,nc);
    if ResampleData == 1
        for k = 1:nVar1
            tempData = [repmat(Data(1,k),64,1); Data(:,k); repmat(Data(end,k),64,1)];   % remove edge effects
            ResampledData(:,k) = resample(tempData, 1, 16);
%      eval([char(VarNames1(k)) '= ResampledData(:,k);']);                                                    % using eval for assigning variable names to the data, be careful
            tempData = [];
        end
    end
%     [nr,~] = size(ResampledData);
    LoopData(:,2:totalVar+1) = ResampledData(5:end-4,1:end);    % row -> (2:end-1) removing edge effects, col -> (2:end) - ignoring timestamp
    clearvars Data ResampledData tempData

    %% Creating a default Variable to inlcude the column headers as tested manually
    % Extract the date string for verification
    folderstr2 = strsplit(string(files2(ii)),'\');
    filestr2 = strsplit(folderstr2(end), 'Hz-');
    datestring2 = filestr2(end);
    if ~isequal(datestring1,datestring2) || ~isequal(numel(files1), numel(files2))
       error('The datetime or the number of files do not match between turbine');
    end

    %% Data import from each file in the directory
%    nfiles = numel(files2);                                                                          % no. of files to be imported
  %  nfiles = 50;                                                                                    % manual entry
    nVar2 = length(VarNames2);                                                                        % no. of variables to be imported
    totalVar = totalVar + nVar2;
    [fid, ~] = fopen(files2{ii},'r');
    NewData = textscan(fid, FormatSpec2, 'Delimiter', ';', 'EmptyValue', NaN,...
               'headerlines', StartRow-1, 'returnonError',0, 'EndOfLine', '\r\n');                 % import data
    fclose(fid);
    cellsize = cellfun(@numel, NewData(1:end-1));
    if ~all(cellsize==mode(cellsize));
        idx = find(cellsize(:)~=mode(cellsize));
        for jc = 1:numel(idx)
            temp = nan(mode(cellsize),1);
            temp(1:cellsize(idx(jc))) = NewData{idx(jc)};
            NewData{idx(jc)} = temp;
            clearvars temp
        end
    end

%     ts = datenum([NewData{:,1}]);                                                                  % time stamp for wind speed
    Data = cell2mat(NewData(1:end-1));                                                        % collecting in one array for writing purposes
    clearvars NewData

    % Check for NaNs and resample the columns
    if (sum(sum(isnan(Data))) > 0) && isempty(find(sum(isnan(Data)) > 0.9*numel(t64)))
      col_nan = find(sum(isnan(Data)) > 0);
      for ic = 1:numel(col_nan)
        method='linear';
        idx_nan = isnan(Data(:,col_nan(ic)));
        if idx_nan(1)==1
            Data(1,col_nan(ic)) = nanmean(Data(1:128,col_nan(ic)));
        elseif idx_nan(end)==1
            Data(end,col_nan(ic)) = nanmean(Data(end-128:end,col_nan(ic)));
        end
        idx_nan = isnan(Data(:,col_nan(ic)));
        Data(idx_nan,col_nan(ic)) = interp1(t64(~idx_nan), Data(~idx_nan,col_nan(ic)), t64(idx_nan), method, 'extrap');
      end % ic #replace nans condn
    elseif ~isempty(find(sum(isnan(Data)) > 0.9*numel(t64)))
        disp('Too many Nans in additional signals XD115');
        nanfiles2 = [nanfiles2; files2(ii)];
    end
    [nr,nc] = size(Data);
    ResampledData = nan(round(nr/16)+8,nVar2);
    if ResampleData == 1
        for k = 1:nVar2
            tempData = [repmat(Data(1,k),64,1); Data(:,k); repmat(Data(end,k),64,1)];   % remove edge effects
            ResampledData(:,k) = resample(tempData, 1, 16);
%            eval([char(VarNames2(k)) '= ResampledData(:,k);']);                                                    % using eval for assigning variable names to the data, be careful
             tempData = [];
        end
    end % 
%     [nr,~] = size(Data);
    LoopData(:,2:totalVar+1) = [LoopData(:,2:nVar1+1), ResampledData(5:end-4,1:end)];
    clearvars Data ResampledData

    %% Creating a default Variable to inlcude the column headers as tested manually
    % Extract the date string for verification
    folderstr3 = strsplit(string(files3(ii)),'\');
    filestr3 = strsplit(folderstr3(end), 'Hz.');
    datestring3 = filestr3(end);
    if ~isequal(datestring1,datestring2, datestring3) || ~isequal(numel(files1), numel(files2), numel(files3))
       error('The datetime or the number of files do not match');
    end

    %% Data import from each file in the directory
    NewData = [];
    Data = [];
    nfiles = numel(files3);                                                                          % no. of files to be imported
    nVar3 = length(VarNames3);                                                                        % no. of variables to be imported
    totalVar = totalVar + nVar3;
    [fid, ~] = fopen(files3{ii},'r');
    NewData = textscan(fid, FormatSpec3, 'Delimiter', ';', 'EmptyValue', NaN,...
               'headerlines', StartRow-1, 'returnonError',0, 'EndOfLine', '\r\n');                 % import data
    fclose(fid);
    ts = datenum([NewData{:,1}]);                                                                  % time stamp for wind speed
    Data = [ts cell2mat(NewData(2:end-1))];                                                        % collecting in one array for writing purposes
    clearvars NewData
%     for k = 1:nVar3
%       eval([char(VarNames3(k)) '=Data(:,k);']);                                                    % using eval for assigning variable names to the data, be careful
%     end
%    [nr,~] = size(Data);
    LoopData(:,1:end) = [Data(:,1), LoopData(:,2:(nVar1+nVar2+1)), Data(:,2:end)];
    clearvars Data

    % Save the workspace variables into a mat file
    if Append == 1
        if ii==1
            save('d:\Main\MATLAB\Current\Wt10_TurbineOperation_full_Torque2018.mat', 'LoopData', '-v7.3');
            count_Append = count_Append + 1;
        else
            mf = matfile('d:\Main\MATLAB\Current\Wt10_TurbineOperation_full_Torque2018.mat', 'Writable', true);
            mf.LoopData(end+1:end+2400, :) = LoopData;
            count_Append = count_Append + 1;
        end % 
    clearvars LoopData
    end % append condn
    fprintf('[%s] - Files processed:%2.0f/%.0f \n', datestr(now,'HH:MM:SS'),ii,nfiles)
end % ii nfiles condn
toc
