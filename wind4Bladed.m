% Script to genererate a Bladed file from wind data

clear; clc; close all
addpath('d:\Delft\SurfDrive\MATLAB\Functions\');
addpath('d:\Delft\SurfDrive\MATLAB\WindModel\Windsim\');

%% Provide default data
h = (185:-15:50)';                                                                                 % Lidar ranges, other related terms Lrange, Le,
Lab = sprintf('Ws%d;',h);                                                                          % wind speed labels for naming
LabN = strsplit(Lab,';');                                                                          % splitting the labels
LabN = LabN(1:10);                                                                                 % create an array of labels
VarWs = sprintf('Wt10_AvPr_D%d_Speed;', h);                                                        % create a string of speed variable names
VarN = strsplit(VarWs,';');                                                                        % split the string names
VarN = VarN(1:10)';                                                                                % create a string array of variable names
VarRws = sprintf('Wt10_AvPr_D%d_RWS0;',h);                                                         % create a string of radial wind speed names
VarR = strsplit(VarRws,';');                                                                       % split the string
VarR = VarR(1:10);                                                                                 % create the rws variable names
time1 = 'Timestamp';                                                                               % time stamp variable name
datum = datenum('2013-12-22 14.00.00','yyyy-mm-dd HH.MM.SS');
[year,month, day, hour, minute, second] = datevec(datum);
DefinitePath = 1;                                                                                  % Definite path as selected using datum, dateExt amd DeffPath var
SmoothData = 1;                                                                                    % option for smoothen the data, '1'-yes, '0'-no
  ResampleData = 0;                                                                                % option for resampling data to remove certain freq,'1'-yes, '0'-no
  winsize = 51;                                                                                   % window size for smoothening
  coeffMa = ones(1,winsize)/winsize;                                                               % coefficiencts for smoothening, without weights
  polyOrd = 1;                                                                                     % polynomial order for smoothening, options-'1'/'2'/'3', pref. lower
  fdelay = (winsize-1)/2;                                                                          % delay due to smoothening
  padvec = ones(fdelay,1);                                                                         % padding vector
    winsize_sp = 201;                                                                              % window size
    coeffMa_sp = ones(1, winsize_sp)/winsize_sp;                                                   % weights
    fdelay_sp = (winsize_sp-1)/2;                                                                  % delay due to smoothening
    padvec_sp = ones(fdelay_sp,1);                                                                 % padding vector
       coeffMa_win = ones(1,11)./11;
       fdelay_win = (length(coeffMa_win)-1)/2;
timevec = (0.25:0.25:600)';                                                                        % time vec for wind speed data
PreFiltering = 1;                                                                                  % "1" prefiltering 'on', "0" prefiltering 'off'
fprintf('[%s]: Default data provided \n', datestr(now,'HH:MM:SS'));

% Import the wind data
String = '^Wt\w*(Speed\w*$|RWS\w*$|RWS_Status$|_Direction$)';                                                    % string to search in the data file, include string to import
StartRow = 4;                                                                                      % start row for data, look in to csv file
NCols = 287;                                                                                      % no. of columns in the data file, look into csv file
dateExt = [int2str(year), '-',num2str(month,'%02d'),'-',num2str(day,'%02d'), ' ',num2str(hour, '%02d'),'.',num2str(minute,'%02d'),'.', '00'];
DeffPath = ['d:\Delft\SurfDrive\Data\data4Hz_WI_AvPr\' num2str(year,'%02d'), '\month ' num2str(month,'%02d'),...
               '\data4Hz.', dateExt, '.csv'] ;                                                     % directory path and file name
if DefinitePath == 0
  % adding the paths for the files to the matrix with name 'files'
  fpath = 'd:\Delft\SurfDrive\Data\data4Hz_WI_AvPr';                                                            % directory path
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
end;

  % creating a default variable to include the column headers as tested manually
  FileId = fopen(DeffPath,'r');
  DefVarNames1 = textscan(FileId, '%s', NCols, 'delimiter', ';', 'HeaderLines',1);                 % all the variables in the csv file
  fclose(FileId);

  % find the indexes of the colum to be imported
  Var = [DefVarNames1{:}];
  ind = regexp(Var, String, 'match');                                                              % regexp to match the variable names to the string input
  indf = find(not(cellfun('isempty', ind)));                                                       % find matching instances
  index = [1;indf];                                                                                % indexes to be importedd
  VarNames1 = Var(index);                                                                          % variables that will be imported
  indfs = indf-1;                                                                                  %
  fmt1 = {'%s'};                                                                                   %
  fmt2 = repmat({'%*f'},1, NCols-1);                                                               %
  fmt2(indfs) = {'%f'};                                                                            %
  FormatSpec1 = [fmt1{:} fmt2{:} '%[^\r\n]'];                                                      % format specification for the import variables
  clearvars ind indf index indfs fmt1 fmt2 FileId NCols sub fls files i j

  %% Data import from each file in the directory
  count = 0;
  Teff_series = [];
  NewData = [];
  Data = [];
  nfiles = numel(files1);                                                                          % no. of files to be imported
%  nfiles = 50;                                                                                    % manual entry
  nVar = length(VarNames1);                                                                        % no. of variables to be imported
%  for ii = 1:nfiles % no. of files to be imported
ii = 1;
    [fid, message] = fopen(files1{ii},'r');
    NewData = textscan(fid, FormatSpec1, 'Delimiter', ';', 'EmptyValue', NaN,...
               'headerlines', StartRow-1, 'returnonError',0, 'EndOfLine', '\r\n');                 % import data
    fclose(fid);
    ts = datenum([NewData{:,1}]);                                                                  % time stamp for wind speed
    Data = [ts cell2mat(NewData(2:end-1))];                                                        % colelcting in one array for writing purposes
    clearvars NewData
    for k = 1:nVar
      eval([char(VarNames1(k)) '=Data(:,k);']);                                                    % using eval for assigning variable names to the data, be careful
    end
    [nr,~] = size(Data);                                                                           % no. of rows in the data variables
fprintf('[%s]: Import data successful \n', datestr(now,'HH:MM:SS'));

    % Prefiltering
if PreFiltering == 1
    fname = files1{ii};                                                                            % file name
    idx_Speed = find(not(cellfun('isempty', regexp(VarNames1,'Speed$','match'))));                 % index for matching the imported variables
    [s,c] = sort_nat(VarNames1(idx_Speed),'descend');                                              % sort the variables in the descending distance acc. to Lrange
    sortidx_Speed = idx_Speed(c);                                                                  % sorted wind speed variables
    nc = length(idx_Speed);                                                                        % no.of columns to be imported
    Ws = nan(nr,nc);                                                                               % works as long as nfiles = 1
    pk = nan(10,nc);                                                                               % "", resets for every file, therefore re-allocate if necessary
    loc = pk; wid = pk; prom = pk;                                                                 % ""
    Ws_raw = zeros(nr, nc);                                                                        % initiation of multiple variables
    wind.mWs = zeros(1,nc);                                                                        %
    wind.stdWs = zeros(1,nc);                                                                      %
    wind.dWs = zeros(nr,nc);                                                                       %
    pad_Ws = zeros(nr+2*length(padvec_sp), nc);                                                    %
    tempWs = zeros(nr+2*length(padvec_sp), nc);                                                    %
    wind.filWs = zeros(nr, nc);                                                                    %
    wind.nMaxima = zeros(1, nc);                                                                   %
    wind.nMinima = zeros(1,nc);                                                                    %
    pkin_max = pk; pkin_min = pk;                                                                  %
    locin_max = loc; locin_min = loc;                                                              %
    widin_max = wid; widin_min = wid;                                                              %
    promin_max = prom; promin_min = prom;                                                          %
    for ic = 1:nc % no. of range gates
      Ws_raw(:,ic) = Data(:,sortidx_Speed(ic));                                                    % raw wind speed measured
      wind.mWs(:,ic) = nanmean(Ws_raw(:,ic));                                                      % mean wind speed in the measrued data
      wind.stdWs(:,ic) = nanstd(Ws_raw(:,ic));                                                     % std deviations
      Ws(:,ic) = FnWsRange(Ws_raw(:,ic), 0.5,2*wind.mWs(:,ic),0);                                  % filter below min & max wind speed values FnWsRange(Ws,WsMin,WsMax)
      idx_nan = (isnan(Ws(:,ic)));                                                                 % find the nan positions
      if idx_nan(1)==1                                                                             % endpoints of time series cannon be interpolated, so
         Ws(1,ic) = nanmean(Ws(1:100,ic));                                                         % assigning the mean over 100 values to the end points
      end                                                                                          %
      if idx_nan(end)==1                                                                           % endpoint cannot be interpolated, correction loop
         Ws(end,ic) = nanmean(Ws(end-100:end,ic));                                                 % assigning the mean over last 100 values
      end                                                                                          %
      idx_nan = (isnan(Ws(:,ic)));                                                                 % find the nan positions
      Ws(idx_nan,ic) = interp1(timevec(~idx_nan),Ws(~idx_nan,ic),timevec(idx_nan), 'nearest');     % interploation using nearest neighbours for nan values
      idx_nan = find(isnan(Ws(:,ic)));                                                             % find nan values again
%        while(~isempty(idx_nan))                                                                  % while all nans are removed
%          Ws(idx_nan,ic) = Ws(idx_nan-1,ic);                                                      % use previous values to fill nans, alternatively use interpolate
%          idx_nan        = find(isnan(Ws(:,ic)));                                                 %
%        end                                                                                       %
      wind.dWs(:,ic) = Ws(:,ic) - wind.mWs(:,ic);                                                  % flucutuations
      pad_Ws(:,ic) = [wind.mWs(:,ic).*padvec_sp; Ws(:,ic); wind.mWs(:,ic).*padvec_sp];             % padding the time series with the mean 0
      tempWs(:,ic) = filter(coeffMa_sp, polyOrd, pad_Ws(:,ic));                                    % smoothened wind speed time series
      wind.filWs(:,ic) = tempWs(2*fdelay_sp+1:end,ic);                                             % truncate the filtered signal
      % gust detection
      [ymax1,~,ymin1,~] = extrema(wind.filWs(:,ic));                                               % extreme points in the time series, all
      wind.nMaxima(:,ic) = length(ymax1);                                                          % no. of maximas
      wind.nMinima(:,ic) = length(ymin1);                                                          % no. of minimas
      [pks_max,locs_max,wids_max,proms_max] = findpeaks(wind.filWs(:,ic),'MinPeakDistance', 30,... %
         'MinPeakProminence',0.5*wind.stdWs(:,ic),'MinPeakWidth',75, 'Annotate', 'extents');       % significant maximas for gust detections
      [pks_min,locs_min,wids_min,proms_min] = findpeaks(-wind.filWs(:,ic),'MinPeakDistance',30,... %
         'MinPeakProminence',0.5*wind.stdWs(:,ic),'MinPeakWidth', 75, 'Annotate', 'extents');      % significant minimas
      pkin_max(1:size(pks_max,1),ic) = pks_max;                                                    % maxima values
      locin_max(1:size(pks_max,1),ic) = locs_max;                                                  % maxima locations
      widin_max(1:size(pks_max,1),ic) = wids_max;                                                  % maxima widths
      promin_max(1:size(pks_max,1),ic) = proms_max;                                                % maxima prominence
      pkin_min(1:size(pks_min,1),ic) = pks_min;                                                    % minima values
      locin_min(1:size(pks_min,1),ic) = locs_min;                                                  % minima locations
      widin_min(1:size(pks_min,1),ic) = wids_min;                                                  % minima widths
      promin_min(1:size(pks_min,1),ic) = proms_min;                                                % minima prominence
      % step detection
      dgx = 0.5;                                                                                   % bin multiplier
      minbin = round((1/dgx)*min(Ws(:,ic)))/(1/dgx);                                               % min bin edge
      maxbin = round((1/dgx)*max(Ws(:,ic)))/(1/dgx);                                               % max bin edge
      gx = minbin-dgx:dgx:maxbin+dgx;                                                              % bin edges
      [Nbin,edges, bin] = histcounts(Ws(:,ic),gx);                                                 % histogram, no. of occurrences in each bin
      wind.stepWs(:,ic) = (edges(bin))';                                                           % converting the time series according to bin values
      % constant step detection                                                                    %
      win = 300;                                                                                   % window size for averaging
      mwinWs = nanmean(reshape(Ws(:,ic),win,[]));                                                  % mean in the window
      stepwinWs = mwinWs(ones(win,1),:);                                                           % step, first 300 values in each column
      stepwinWs = stepwinWs(:)';                                                                   % reshape concatenate
      wind.stepWs2(:,ic) = stepwinWs;                                                              % concatenated wind speed
      % Another method for mean steps detection for dynamic inflow purposes
      [ipt,~] = findchangepts(Ws(:,ic),'Statistic','mean','MinThreshold',200,'MinDistance',480);   % find abrupt changes in wind speed
      wind.nSteps(ic,ii) = numel(ipt);                                                             % no. of steps in wind speed
    end % ic, range gates condn
    wind.nSteps_ui(ii) = mode(wind.nSteps(:,ii));                                                  % no. of step changes in wind speed in a file
    clearvars idx_Speed idx_nan ymax1 imax ic ipt stepwinWs pks_max pks_min locs_min locs_max wids_min wids_max proms_min proms_max ymax1 ymin1 Ws_raw...
         win mwinWs pk loc wid prom count Data pad_Ws tempWs
fprintf('[%s]: pre-processing data successful \n', datestr(now,'HH:MM:SS'));
end


% Add the parameters required for writing to bladed
fprintf('[%s]: Adding Lidar data to turbulent wind field \n', datestr(now,'HH:MM:SS'));
filename = 'data64Hz.2013-12-22 14.00.00';
T = 600; % Time series, total time [s]
dt = 0.25; % time step [s]
t = (dt:dt:T)';
Ly = 200;  % domain width [m]
Lz = 200;  % domain height [m]
Ny = 7;   % grid points in y-direction/lateral [-]
Nz = 7;   % grid points in z-direction/vertical [-]
Uhub = nanmean(Ws(:,end)); % hub height wind speed [m/s]
% Uhub = 10; % hub height averaged wind speed [m/s]
H = 100; % hub height [m]
dx = dt*Uhub; % distance between two grid points [m]
Nfft = T/dt; % no. of fft terms
x = t.*Uhub; % length of the domain in x-direction [m]
Blindzone = 50; % blind zone of the lidar in front of the wt [m]
bufferlength = 150;   % bufferlength in case of negative distances [m]
nBeams = 5; % no., of Lidar beams
Lx = x(end)+Blindzone + bufferlength;
theta = 15; % cone angle of Lidar [deg]
rangeLidar = (185:-15:Blindzone); % lidar range gates
rangeFull = (Lx:-dx:dx); % complete domain range
rangeBeam = (Lx-dx:-dx:Blindzone); % range of Lidar beams
range = [ repmat(rangeLidar,1,4), rangeLidar(1:8)]; %
range = rangeLidar;
r = Lx - rangeLidar;
dy = Ly/(Ny-1);
dz = Lz/(Nz-1);
df = 1/T;
fn = df/2;
Fs = Nfft/T;
f = (df:df:0.5*Fs);
kw = 2*pi*f/Uhub;

%% Domain
y = linspace(Ly/2,-Ly/2, Ny);
z = H + linspace(-Lz/2, Lz/2, Nz);
[X,Y,Z] = ndgrid(x,y,z);
pos = [X(:) Y(:) Z(:)];                                                   % concatenated matrix
DT  = delaunay(pos);                                                      % delaunay matrix for relative distances
fprintf('[%s]: Domain and relative positions created \n', datestr(now,'HH:MM:SS'));

% Assign grid points to Lidar measurements
nt = numel(t);
t0 = 1:5:nt; % timestamp of beam0
t1 = 2:5:nt; % timestamp of beam1
t2 = 3:5:nt; % timestamp of beam2
t3 = 4:5:nt; % timestamp of beam3
t4 = 5:5:nt; % timestamp of beam4

WsMin = 2;
WsMax = 25;
rws50 = FnWsRange(Wt10_AvPr_D50_RWS, WsMin, WsMax,0);
rws65 = FnWsRange(Wt10_AvPr_D65_RWS, WsMin, WsMax,0);
rws80 = FnWsRange(Wt10_AvPr_D80_RWS, WsMin, WsMax,0);
rws95 = FnWsRange(Wt10_AvPr_D95_RWS, WsMin, WsMax,0);
rws110 = FnWsRange(Wt10_AvPr_D110_RWS, WsMin, WsMax,0);
rws125 = FnWsRange(Wt10_AvPr_D125_RWS, WsMin, WsMax,0);
rws140 = FnWsRange(Wt10_AvPr_D140_RWS, WsMin, WsMax,0);
rws155 = FnWsRange(Wt10_AvPr_D155_RWS, WsMin, WsMax,0);
rws170 = FnWsRange(Wt10_AvPr_D170_RWS, WsMin, WsMax,0);
rws185 = FnWsRange(Wt10_AvPr_D185_RWS, WsMin, WsMax,0);

beam0 = [rws50(t0),rws65(t0),rws80(t0),rws95(t0),rws110(t0),rws125(t0),rws140(t0),rws155(t0), rws170(t0),rws185(t0)];
beam1 = [rws50(t1),rws65(t1),rws80(t1),rws95(t1),rws110(t1),rws125(t1),rws140(t1),rws155(t1),rws170(t1),rws185(t1)];
beam2 = [rws50(t2),rws65(t2),rws80(t2),rws95(t2),rws110(t2),rws125(t2),rws140(t2),rws155(t2),rws170(t2),rws185(t2)];
beam3 = [rws50(t3),rws65(t3),rws80(t3),rws95(t3),rws110(t3),rws125(t3),rws140(t3),rws155(t3),rws170(t3),rws185(t3)];
beam4 = [rws50(t4),rws65(t4),rws80(t4),rws95(t4),rws110(t4),rws125(t4),rws140(t4),rws155(t4),rws170(t4),rws185(t4)];

dbeam0 = [Wt10_AvPr_D50_DRWS(t0),Wt10_AvPr_D65_DRWS(t0),Wt10_AvPr_D80_DRWS(t0),Wt10_AvPr_D95_DRWS(t0),Wt10_AvPr_D110_DRWS(t0),Wt10_AvPr_D125_DRWS(t0),Wt10_AvPr_D140_DRWS(t0),rws155(t0), rws170(t0),rws185(t0)];
dbeam1 = [Wt10_AvPr_D50_DRWS(t1),Wt10_AvPr_D65_DRWS(t1),Wt10_AvPr_D80_DRWS(t1),Wt10_AvPr_D95_DRWS(t1),Wt10_AvPr_D110_DRWS(t1),Wt10_AvPr_D125_DRWS(t1),Wt10_AvPr_D140_DRWS(t1),rws155(t1),rws170(t1),rws185(t1)];
dbeam2 = [Wt10_AvPr_D50_DRWS(t2),Wt10_AvPr_D65_DRWS(t2),Wt10_AvPr_D80_DRWS(t2),Wt10_AvPr_D95_DRWS(t2),Wt10_AvPr_D110_DRWS(t2),Wt10_AvPr_D125_DRWS(t2),Wt10_AvPr_D140_DRWS(t2),rws155(t2),rws170(t2),rws185(t2)];
dbeam3 = [Wt10_AvPr_D50_DRWS(t3),Wt10_AvPr_D65_DRWS(t3),Wt10_AvPr_D80_DRWS(t3),Wt10_AvPr_D95_DRWS(t3),Wt10_AvPr_D110_DRWS(t3),Wt10_AvPr_D125_DRWS(t3),Wt10_AvPr_D140_DRWS(t3),rws155(t3),rws170(t3),rws185(t3)];
dbeam4 = [Wt10_AvPr_D50_DRWS(t4),Wt10_AvPr_D65_DRWS(t4),Wt10_AvPr_D80_DRWS(t4),Wt10_AvPr_D95_DRWS(t4),Wt10_AvPr_D110_DRWS(t4),Wt10_AvPr_D125_DRWS(t4),Wt10_AvPr_D140_DRWS(t4),rws155(t4),rws170(t4),rws185(t4)];

fprintf('[%s]: Beam measurements categorized \n', datestr(now,'HH:MM:SS'));

% beam positions in general
U0 = nan(Nfft*Ny*Nz,1);
U1=U0; U2 = U0; U3 = U0; U4 = U0; W1 = U0; W3 = U0; V2 = U0; V4 = U0;

for it = 1:Nfft/nBeams
  if it==1
    b0 = [r; +0*range;                  H + 0*range];
    b1 = [r; +0*range;                  H + abs(range)*tand(theta) ];
    b2 = [r; -abs(range)*tand(theta);   H + 0*range ];
    b3 = [r; +0*range;                  H-abs(range)*tand(theta)];
    b4 = [r; +abs(range)*tand(theta);   H + 0*range ];
    fprintf('[%s]: Beam positions defined \n', datestr(now,'HH:MM:SS'));
  else
    b0(1,:) = b0(1,:) - nBeams*dx;
    b1(1,:) = b1(1,:) - nBeams*dx;
    b2(1,:) = b2(1,:) - nBeams*dx;
    b3(1,:) = b3(1,:) - nBeams*dx;
    b4(1,:) = b4(1,:) - nBeams*dx;
  end
  [k0(:,it), ~] = dsearchn(pos,DT, b0');                                        % nearest point search in the grid for Lidar measured ranges
  [k1(:,it), ~] = dsearchn(pos,DT, b1');                                        % nearest point search
  [k2(:,it), ~] = dsearchn(pos,DT, b2');                                        % nearest point search
  [k3(:,it), ~] = dsearchn(pos,DT, b3');                                        % nearest point search
  [k4(:,it), ~] = dsearchn(pos,DT, b4');                                        % nearest point search

  % Assigning the wind speeds to the wind field
  U0(k0(:,it)) = (beam0(it,:)./cosd(0));                                       % re-fitting Lidar measurements in Longitudinal direction assuming frozen turbulence
  U1(k1(:,it)) = (beam1(it,:)./cosd(theta));
  U2(k2(:,it)) = (beam2(it,:)./cosd(theta));
  U3(k3(:,it)) = (beam3(it,:)./cosd(theta));
  U4(k4(:,it)) = (beam4(it,:)./cosd(theta));
end
fprintf('[%s]: Added Lidar measurements to nearest grid positions \n', datestr(now,'HH:MM:SS'));

% Reshaping into the domain size
% mean of respective beams
mU0 = nanmean(U0);
mU1 = nanmean(U1);
mU2 = nanmean(U2);
mU3 = nanmean(U3);
mU4 = nanmean(U4);
% standard deviation in the directions
uu0 = nanstd(U0);
uu1 = nanstd(U1);
uu2 = nanstd(U2);
uu3 = nanstd(U3);
uu4 = nanstd(U4);
% fluctuations in the directions
u0 = U0 - mU0;
u1 = U1 - mU1;
u2 = U2 - mU2;
u3 = U3 - mU3;
u4 = U4 - mU4;

U = nan(Nfft*Ny*Nz,1);    % mean wind field
u = nan(Nfft*Ny*Nz,1); % longitudinal component
v = nan(Nfft*Ny*Nz,1); % lateral component
w = nan(Nfft*Ny*Nz,1); % vertical component
uu  = nan(Nfft*Ny*Nz,1); % longitudinal fluctuations component
vv  = nan(Nfft*Ny*Nz,1); % lateral fluctuations component
ww  = nan(Nfft*Ny*Nz,1); % vertical fluctuations component

u(~isnan(U0)) = U0(~isnan(U0));
u(~isnan(U1)) = U1(~isnan(U1));
u(~isnan(U2)) = U2(~isnan(U2));
u(~isnan(U3)) = U3(~isnan(U3));
u(~isnan(U4)) = U4(~isnan(U4));
U = sqrt(u.^2 + v.^2 +w.^2);

uu(~isnan(u0)) = u0(~isnan(u0));
uu(~isnan(u1)) = u1(~isnan(u1));
uu(~isnan(u2)) = u2(~isnan(u2));
uu(~isnan(u3)) = u3(~isnan(u3));
uu(~isnan(u4)) = u4(~isnan(u4));
vv(~isnan(v2)) = v2(~isnan(v2));
vv(~isnan(v4)) = v4(~isnan(v4));
ww(~isnan(w1)) = w1(~isnan(w1));
ww(~isnan(w3)) = w3(~isnan(w3));

U0 = reshape(U0 ,Nfft, Ny, Nz);
U1 = reshape(U1, Nfft,Ny, Nz);
U2 = reshape(U2, Nfft,Ny, Nz);
U3 = reshape(U3, Nfft,Ny, Nz);
U4 = reshape(U4, Nfft,Ny, Nz);
u = reshape(u, Nfft, Ny, Nz);
v = reshape(v, Nfft, Ny, Nz);
w = reshape(w, Nfft, Ny, Nz);
uu = reshape(uu, Nfft, Ny, Nz);
vv = reshape(vv, Nfft, Ny, Nz);
ww = reshape(ww, Nfft, Ny, Nz);

fprintf('[%s]: Lidar measurements integrated into the domain \n', datestr(now,'HH:MM:SS'));
error('manual stop');
% longitudinal variations only
u = zeros(size(X));
v = u;
w = u;
u(1:numel(x),6,6) = Ws(:,end);

writebladed('data64Hz.2013-12-22 14.00.00.wnd',u,v,w,x,y,z,U)

