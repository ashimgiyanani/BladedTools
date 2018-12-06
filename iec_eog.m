% Script to create a extreme operating gust in Matlab
% References: IEC2010 IEC 61400-1 Ed.3 Wind turbines - Design requirements
% Terms used as used in the IEC report

% Initializations
write_wnd = 'writeBLgrid'; % options 'writebladed' or 'writeBLgrid'
validate_wnd = 1;
zhub = 126;
edition = 2;

%% Values according to Table 1 in IEC-61400-1 Ed.3
if edition == 3
  WTC_param.classes = [1, 2, 3]; % wind turbine classes for conditions excl. offshore, and complex sites
  WTC_param.Vref = [50, 42.5, 37.5]; % reference wind speed average over 10 min [m/s]
  WTC_param.Iref = [0.16, 0.14, 0.12]; % expected value of the turbulence intensity at 15 m/s [fraction]
  WTC_param.category = ['A', 'B', 'C']; % A - high turbulence, B - medium turbulence, C -  low turbulence
elseif edition == 2
  WTC_param.classes = [1, 2, 3, 4]; % wind turbine classes for conditions excl. offshore, and complex sites
  WTC_param.Vref = [50, 42.5, 37.5, 30]; % reference wind speed average over 10 min [m/s]
  WTC_param.Vave = [10, 8.5, 7.5, 6]; % reference wind speed average over 10 min [m/s]
  WTC_param.Iref = [0.18,0.16]; % expected value of the turbulence intensity at 15 m/s [fraction]
  WTC_param.a = [2,3];
  WTC_param.category = ['A', 'B']; % A - high turbulence, B - medium turbulence, C -  low turbulence
  N = 50;
end

D = 141; % rotor diameter [m]
WTC = 1; % wind turbine class
turbCat = 'A'; %  
idx_wtc = find(WTC_param.classes==WTC); % index to choose other IEC parameters
idx_turbCat = find(WTC_param.category == turbCat); % index to choose turbulence category
Iref = WTC_param.Iref(idx_turbCat);  % reference turbulence intensity, based on wtc
Vref = WTC_param.Vref(idx_wtc); % reference wind speed based on wtc
a = WTC_param.a(idx_turbCat);  % reference a, based on wtc



%% Extreme wind speed model (EWM)
T = 10.5;
dt = 0.1;
t = (dt:dt:T)';
Vhub = 25.*ones(size(t)); % hub height wind speed [m/s], e.g.
% In the steady extreme wind model, allowance for short-term deviations from the mean wind direction shall be made by assuming constant yaw misalignment in the range of ±15º.
zhub = 30; % hub height [m] e.g.
z = 30; 
Ve50 = 1.4*Vref*(z/zhub)^0.11
Ve1 = 0.8*Ve50;
% Turbulent extreme wind model
sig1 = 0.11*Vhub



%% Extreme operating wind gust
D = 141;
dt = 0.01;
Tfinal = 80.36;
t = (dt:dt:Tfinal)';
Vhub = 8.3.*ones(size(t)); % hub height wind speed [m/s], e.g.
zhub = 126; % hub height [m], e.g.
z = (linspace((zhub - D/2), (zhub + D/2), 21));

if edition == 3
    % Time of the event
    T = 10.5;
    % Calculating Lambda1
    if zhub <= 60
      lambda1 = 0.7.*z; % turbulence scale parameter [m]
    elseif zhub > 60
      lambda1 = 42; % turbulence scale parameter [m]
    end
    % Calculating sigma1
    b = 5.6;  
    sig1 = Iref.*(0.75.*Vhub + b); % for NTM, turbulence standard deviation given by 90% quantile of Vhub
    % Calculating Vgust
    Ve50 = 1.4.*Vref.*((z./zhub).^0.11); % extreme wind speed with a recurrence period of 50 years
    Ve1 = 0.8.*Ve50; % extreme wind speed with a recurrence period of 50 years
    Vgust1 = 1.35.*(Ve1 - Vhub); % 
    Vgust2 = 3.3.*(sig1./(1+0.1.*(D./lambda1)));
    Vgust = min(Vgust1, Vgust2);

elseif edition == 2
   % Calculating Lambda1
    if zhub <= 30
      lambda1 = 0.7.*z; % turbulence scale parameter [m]
    elseif zhub > 30
      lambda1 = 21; % turbulence scale parameter [m]
    end
    % Calculating sigma1
    sig1 = Iref.*(15 + a.*Vhub)./(a+1); % for NTM, turbulence standard deviation given by 90% quantile of Vhub
    if N==1
        betaN= 4.8; % factor for 1 year EOG N=1
        T = 10.5;   % time of the event
    elseif N==50
        betaN= 6.4; % factor for 1 year EOG N=50
        T = 14;     % time of the event
    end
    Vgust = betaN.*(sig1./(1+0.1.*(D./lambda1)));
end

% Calculate the wind speed acc. to wind profile, same for ed.2 and ed.3
Alpha = 0.2; % power law exponent
Vz = Vhub.*((z./zhub).^Alpha); % wind profile, wind speed as a function of z 

% Simulation event details
stepStart = 55; % start of the event in time series
stepEnd = stepStart + T; % end of the event in the timeseries
idx_Tstart = find(t == stepStart+dt); % index of the event starting time
idx_Tend = find(t == stepEnd); % index of the event ending time
Vzt = Vz; 
Vzt(idx_Tstart:idx_Tend,:) = Vz(idx_Tstart:idx_Tend,:) - 0.37.*Vgust(idx_Tstart:idx_Tend).*sin(3.*pi.*t(idx_Tstart:idx_Tend)./T).*(1 - cos(2.*pi.*t(idx_Tstart:idx_Tend)./T));

figure;
plot(t, Vzt(:,11), 'LineWidth', 2);
grid on;
% axis([t(idx_Tstart) t(idx_Tend) min(Vzt)-2 max(Vzt)+2]);
title('extreme operating gust (EOG)');
xlabel('Time [s]');
ylabel('Wind speed [m/s]');

%%  Extreme wind shear (EWS)
t=(0:dt:T)';
T = 12;
betaN=6.4;
stepStart = 0;
stepEnd = stepStart + T;
Vhub= 8.3.*ones(size(t));
idx_Tstart = find(t == stepStart) + 1; % index of the event starting time
idx_Tend = find(t == stepEnd); % index of the event ending time
z = (linspace((zhub - D/2), (zhub + D/2), 21));
% Calculate the wind speed acc. to wind profile, same for ed.2 and ed.3
Alpha = 0.2; % power law exponent
Vz = Vhub.*((z./zhub).^Alpha); % wind profile, wind speed as a function of z 
% positive EWS
V_ews = (((z - zhub)./D).*(2.5 + 0.2.*betaN.*sig1(idx_Tstart:idx_Tend).*((D./lambda1).^0.25))).*(1 - cos(2.*pi.*t(idx_Tstart:idx_Tend)./T));
wind.ed2.ews.Vzt_p = Vz;
wind.ed2.ews.Vzt_p(idx_Tstart:idx_Tend,:) = Vz(idx_Tstart:idx_Tend,:) + V_ews(1:end,:);
% negative EWS
wind.ed2.ews.Vzt_n = Vz;
wind.ed2.ews.Vzt_n(idx_Tstart:idx_Tend,:) = Vz - V_ews;
plot([wind.ed2.ews.Vzt_p(:,1), wind.ed2.ews.Vzt_p(:,11), wind.ed2.ews.Vzt_p(:,end)])
error('manual stop')
%% Extreme direction change
dt = 0.1; % time step
Tedc = 6; % duration of extreme direction change
t = (-5+dt:dt:15); % total time 
b = 5.6; % 
Vsel = 25; % hub wind speed selected for transient EDC
V = (0:1:40)'; % wind speed steps  
idx_Vsel = find(V == Vsel);
Vhub = V*ones(size(t));
sig1 = Iref.*(0.75.*Vhub + b); % for NTM, turbulence standard deviation given by 90% quantile of Vhub
zhub = 30;
if zhub <= 60
  lambda1 = 0.7*zhub; % turbulence scale parameter [m]
elseif zhub > 60
  lambda1 = 42; % turbulence scale parameter [m]
end 
Theta_e = 4.*atand(sig1./(Vhub.*(1+0.1.*(D./lambda1)))); % extreme direction change magnitude, limited to +-180°
idx_edc_start = find(t==0); % index of start of extreme direction change
idx_edc_end = find(t == Tedc); % index of end of extreme direction change
Thetat = zeros(size(V,1), size(t,2)); % initialization
Thetat(:,idx_edc_start:idx_edc_end) = 0.5.*Theta_e(:,idx_edc_start:idx_edc_end).*(1 - cos(pi.*(ones(size(V))*t(idx_edc_start:idx_edc_end))./Tedc));
Thetat(:,idx_edc_end+1:end) = Theta_e(:,idx_edc_end+1:end);

figure;
plot(V,Theta_e./2, 'LineWidth',2);
hold on;
plot(V, -Theta_e./2, 'lineWidth', 2);
hold off;
grid on;
title('extreme direction change (EDC)');
xlabel('Wind speed  [m/s]');
ylabel('extreme Wind direction change magnitude [deg] ');

figure;
plot(t, Thetat(idx_Vsel,:), 'LineWidth',2);
grid on;
% axis([t(1) t(end) 0 40]);
title('extreme direction change transient');
xlabel('Time [s]');
ylabel('extreme direction change transient [deg]');

%% Extreme operating gust with extreme direction change (ECD)
Vcg = 15; % magnitude of the extreme coherent gust
Tecd = 10; % rise time
dt = 0.1;
Vsel = 25;
zhub = 30; % hub height [m], e.g.
z = 30;
Alpha = 0.2; % power law exponent
t = (-5+dt:dt:15); % total time 
idx_ecd_start = find(t==0);
idx_ecd_end =find(t==Tecd);
V = (0:1:40)'; % wind speed steps  
idx_Vsel = find(V == Vsel);
Vhub = V*ones(size(t));
z = 30;
Vz = Vhub.*((z/zhub)^Alpha); % wind profile, wind speed as a function of z 
Vzt = Vz;
Vzt(:,idx_ecd_start:idx_ecd_end) =  Vz(:,idx_ecd_start:idx_ecd_end) + 0.5.*Vcg.*(1 - cos(pi.*(ones(size(V))*t(idx_ecd_start:idx_ecd_end))./Tecd)); % 
Vzt(:,idx_ecd_end+1:end) = Vz(:,idx_ecd_end+1:end) + Vcg; % 

idx_V4 = find(Vhub < 4); 
idx_Vr = find(Vhub>=4 & Vhub<Vref);
Theta_cg = zeros(size(Vhub));
Theta_cg(idx_V4) = 180;
Theta_cg(idx_Vr) = 720./Vhub(idx_Vr);

Thetat = zeros(size(V,1), size(t,2));
Thetat(:,idx_ecd_start:idx_ecd_end) = 0.5.*Theta_cg(:,idx_ecd_start:idx_ecd_end).*(1-cos(pi.*(ones(size(V))*t(idx_ecd_start:idx_ecd_end))./Tecd));
Thetat(:,idx_ecd_end+1:end) = Theta_cg(:,idx_ecd_end+1:end);
t = (dt:dt:20)
figure();
plot(V, Theta_cg, 'LineWidth', 2);
title('Direction change (ECD)')
xlabel('Wind speed  [m/s]');
ylabel('Wind direction change magnitude [deg]');

figure();
plot(t, Thetat(idx_Vsel,:), 'LineWidth',2)
title('Direction change transient (ECD)')
xlabel('Time  [s]');
ylabel('Wind direction change magnitude [deg]');

figure();
plot(t, Vzt(idx_Vsel,:), 'LineWidth',2)
title('extreme coherent gust  amplitude (ECD)');
xlabel('Time [s]');
ylabel('Wind speed [m/s]');


if strmatch(write_wnd, 'writeBLgrid')
    %% Write to Bladed using writeBLgrid
    FileName = 'TrialBladed_v1_blgrid.wnd';
    Nc = 3; % no. of  components
    Uhub = 25; % hub height mean wind speed
    vzt = Vzt(idx_Vsel,:) - Vsel; % 
    u = vzt.*cosd(Thetat(idx_Vsel,:));
    v = vzt.*sind(Thetat(idx_Vsel,:));
    w = vzt.*0;
    uc = zeros(Nt,Ny,Nz);
    vc = zeros(Nt,Ny,Nz);
    wc = zeros(Nt,Ny,Nz);
    for iy = 1:Ny
        for iz = 1:Nz
          uc(:,iy, iz) = u';
          vc(:,iy, iz) = v';
          wc(:,iy, iz) = w'; 
        end
    end
    uvw = zeros(Nt, Nc, Ny, Nz);
    uvw(:,1,:,:) = uc;
    uvw(:,2,:,:) = vc;
    uvw(:,3,:,:) = wc; 
    yr = (0:30:300)';
    zr = (0:30:300)';
    xr = Vsel.*t;
    zhub = 125.9;
    [valc, idxc] = min(abs(zr-zhub));
    zr(idxc) = zhub
    Hz = find(zr==zhub);
    Hy = find(yr==max(yr)/2);
    [X Y Z] = ndgrid(xr, yr, zr);
    [Nt, Nuvw] = size(uvw)
    Ny = length(yr);
    Nz = length(zr);
    
    % scaling the time series, if sigmau==1
    % fprintf('[%s] - Scaling time series \n', datestr(now,'HH:MM:SS'))
    % u = u * sigmau/std(u(:));
    % v = v * sigmav/std(u(:));
    % w = w * sigmaw/std(u(:));
    % 
    % Include the wind shear
    % fprintf('[%s] - Including wind shear \n', datestr(now,'HH:MM:SS'))
    % Alpha = 0.096*log10(z0) + 0.016*(log10(z0)^2)+0.24;
    % u_shear = Uhub.*((zgv/H).^Alpha);
    
     
%     velocity = zeros(Nt*Nuvw*Ny*Nz,1);
%     velocity = reshape(velocity, Nt, Nuvw, Ny, Nz); 
%     velocity(:,:,Hy, Hz) = uvw;
    % velocity(Nt*Nuvw*Hy*Hz:Nt*Nuvw*(Hy+1)*(Hz+1)) = uvw; % 4D array (time, 3D-wind comp, y, z)
    dy = mean(diff(yr));
    dz = mean(diff(zr));
    dt = 0.1;
    zOffset = 125.9;
    z0 = 0.01;
    Clockwise = 1;
    Ubar = 25;
    I_u = 12;
    I_v = 9.6;
    I_w = 6; 
    SummVars = [zhub; Clockwise; Ubar; I_u; I_v; I_w];
    writeBLgrid(FileName, uvw, dy, dz, dt, zOffset, z0, SummVars);
end

%clearvars -except FileName
% if ValidateData == 1

% end

if strmatch(write_wnd,'writebladed')  
    % use writebladed.m
    FileName = 'TrialBladed_v2_writebl';
    t = (dt:dt:20);
    Ly = 300;
    Lz = 300;
    Ny = 11;
    Nz = 11;
    H = 150;
    Ny = 11;
    Nz = 11;
    yr = linspace(Ly/2,-Ly/2, Ny);
    zr = H + linspace(-Lz/2, Lz/2, Nz);
    Uhub = Vsel;
    xr = t.*Vsel;
    Nt = numel(xr);
    [X Y Z] = ndgrid(xr,yr,zr);
    vzt = Vzt(idx_Vsel,:) - Uhub;
    u = vzt.*cosd(Thetat(idx_Vsel,:));
    v = vzt.*sind(Thetat(idx_Vsel,:));
    w = vzt.*0;
    uc = zeros(Nt,Ny,Nz);
    vc = zeros(Nt,Ny,Nz);
    wc = zeros(Nt,Ny,Nz);
    for iy = 1:Ny
        for iz = 1:Nz
          uc(:,iy, iz) = u';
          vc(:,iy, iz) = v';
          wc(:,iy, iz) = w'; 
        end
    end
    writebladed(FileName,uc,vc,wc,xr',yr',zr',Uhub)
end

if (validate_wnd == 1)
  filename = which(strcat(FileName, '.wnd')); 
  [Folder, baseFileName, extension] = fileparts(filename);
  [u,v,w,x,y,z,U,t,Nx, Ny, Nz, delta_x, delta_y, delta_z, uvw]=readbladed(filename);                    
  spd = U + sqrt(u.^2 + v.^2 + w.^2);
  
  %% Visualizing the 3d wind field
  cla;
  [m,n,p] = size(u);
  [Cx, Cy, Cz] = ndgrid(1:4:m,1:1:n,1:1:p);
  h  = coneplot(u, v, w, Cx, Cy, Cz);
  set(h,'EdgeColor', 'none')
  axis tight equal
  view(37,32)
  box on
  colormap(hsv)
  light

  figure;
  plot(x,spd(:,6,6), 'LineWidth',2);
  xlabel('Length [m]');
  ylabel('Wind speed u component [m/s]');
  
  figure;
  [X,Y, Z] =  meshgrid(xr, yr, zr);
  spd_p = permute(spd,[2 1 3] );
  hiso = patch(isosurface(X,Y,Z, spd_p));
  isonormals(X,Y,Z, spd, hiso);
  hiso.FaceColor = 'red';
  hiso.EdgeColor = 'none';
end
