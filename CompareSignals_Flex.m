% Script to compare signals from Flex output 1 and Flex output 2 i.e. compare Flex outputs with different settings

importFlex = 0;
FlexFile1 = 'c:\Users\ashim.giyanani\OneDrive - Windwise GmbH\Windwise\Flex\ExtremeWindShear_rev04_126mHH\Positive_EWM_Inflow_min5deg.txt';
FlexFile2 = 'c:\Users\ashim.giyanani\OneDrive - Windwise GmbH\Windwise\Flex\ExtremeWindShear_rev04_126mHH\Positive_EWM_DOF_deactivated_Inflow_min5deg.txt';

% import
Delimiter = '\t'; % delimiter in the file to be read
nHead = 14; % number of lines before header names
String = 'VuHc|Teta2|OmRot|P_Gen|FzHR|FzYB|MyTB|MzHR|FzH2|MyYB|MyHR|MxHF|Mx*|My*|Mz*|Fx*|Mfric*|\wH1$|M*\w11$';
fid = fopen(FlexFile1)
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
[fid, message] = fopen(FlexFile1,'r');
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

idx_Time = find(Time == 60.16)
Mx = Mx11(idx_Time) + Fz11(idx_Time).*y_shear_centre; % moment around the shear centre for a radial location or for a particular time
dr = [0;diff(Length)]; 
local_Twist = rad2deg(dr.*Mx.*10^3./Tor_S);
cum_locTwist = cumsum(local_Twist);


