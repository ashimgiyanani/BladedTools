%% Script to create the Synthetic wind file for Feike for wind files coming as close to the measured files
clear; clc; close all;

%% Providing default data
Delimiter = ';';                                                      % delimiter for reading file
StartRow = 2;                                                         % start row for reading the file
File = 'd:\Main\MATLAB\WindModel\ModelData\EvolData_Torque2212130420.txt';%%% Change
SaveFig = 0;                                                          % savefigures yes=1, no= 0
Append = 0;                                                           % append files into the text file, yes(1)/no(0)
fid = fopen(File,'r');
line1 = fgetl(fid);                                                   % get the headers
NCols = length(find(line1 == Delimiter))+1;			                      %  no. of columns based on the delimiter defined
FormatSpec = [repmat('%f', 1, NCols) '\r\n'];                         % Formatspec to read the file
frewind(fid);                                                         % going back to the start of file
Data = textscan(fid, FormatSpec,'Delimiter', Delimiter,...
         'EmptyValue', NaN,'MultipleDelimsAsOne',1,...
           'Headerlines', StartRow-1,'returnOnError', 1);             % reading the data file with config.
fclose(fid);                                                          % close the opend file
Data = [Data{:}];                                                     % converting the array to mat
DefVarNames = strsplit(line1,Delimiter);                              % extracted var names
for k = 1:length(DefVarNames)-1
eval([char(DefVarNames(k)) '= Data(:,k);']);                          % assign variable names to the data
end;
%% Initialization box domain, turbulence field properties and so on
yr = [100];								  																		      % width of the domain
zr = [100];								  																		      % height of the domain
deltaX = 15;										 																			% longitudinal distance between two measnt pts
h = [185:-deltaX:50]';					  																		% Range gates
fmax = 12;                                       											% sampling frequency, why 12 Hz????
deltat = 1/fmax;																											% sampling time
N = (10+1).*60.*fmax;																								  % approximately the length of the time series
Z = 100;																															% Hub height
Zo = 0.5;																															% Surface roughness length
% Extracting the mean and devaition from the Lidar measurements
mWs = [mWs1; mWs2(end)];                                              % mean wind speed Lidar measurements
stdWs = [stdWs1; stdWs2(end)];                                        % std dev Lidar wind speed measurements
StartTime = datetime('2013-12-22 04:30:00.00',...
                   'Format', 'yyyy-MM-dd hh:mm:ss.SS');               % start time of the simulation
Timestamp = [linspace(StartTime, StartTime+minutes(11)...
                                -seconds(deltat), N)]';                  % Creating the time stamp for the time series
numT = datenum(Timestamp);
U = mWs;																															% mean wind speed
fricU = 0.4.*U./(log(Z./Zo));																					% friction velocity
sig = sqrt(0.57.*fricU.^2);																						% wind velocty variation
sigma = sig;
sigma = stdWs;																												% based on ti=sigma./U, IEC
zhub = Z;
tau = h./U;																											      % mean time delay assuming Taylor's hypothesis
tauLag = ceil(tau./deltat);																						% number of time lags for the time delay
Na = [1:tauLag:tauLag.*length(h)]';
Nb = [N:tauLag:N+tauLag.*(length(h)-1)]';
for rr = 1:length(mWs)
%[t,UC, Amp, Phi]=Modwind0(yr,zr,U,sigma,N,deltat,fmax, zhub, Zo);
[t,UC(:,rr), Amp(:,rr), Phi(:,rr)]=wind0(yr,zr,U(rr),...
                                          sigma(rr),N,deltat,fmax);   % Veer's simulation for wind field
UCp(:,rr) = U(rr) + UC(:,rr);                                         % Replication of Lidar ranges measurement based on swf
UCa(1:N-1*60.*fmax,rr) = UCp((tauLag(rr)+1:N-1.*60.*fmax+tauLag(rr))', rr);      % Shifting the time by time lag based on distance h
t = t(1:N-1*60.*fmax);
end;                                                                  % 

% Write  the data to a file
if Append == 1
h = [185:-15:50]';                                                    % Lidar ranges
Lab = sprintf('RWS0_%dm;', h);                                        % Variable name to be assigned
Lab(end) = [];
LabN = strsplit(Lab,';');
LabN = LabN(1:10);                                                    % VAriable names as a matrix
TS = (datestr(numT,'yyyy-mm-dd HH:MM:SS.FFF'));                       % note that the formatOut is diff. for datestr fn.
fid = fopen('Bladed_swfSim2212130420_seedv2.txt','w');                % text file for data
fprintf(fid, ['Timestamp' ';' Lab '\r\n']);                           % writing to text file
fSpec=['%2.2f;%2.2f;%2.2f;%2.2f;%2.2f;%2.2f;%2.2f;'...
                      '%2.2f;%2.2f;%2.2f \r\n'];                      % Formatspec for the writing to file
for i = 1:length(UCa)
fprintf(fid, '%s;', (TS(i,:)));
fprintf(fid, fSpec,  UCa(i,:));                                       % writing the data
end;
fclose(fid);
end;



