function [u,v,w,x,y,z,U,t, Nx, Ny, Nz, delta_x, delta_y, delta_z, uvw]=readbladed(filename)
% syntax function [u,v,w,x,y,z,U,tNx,Ny,Nz,delta_x,delta_y,delta_z,ug,vg,wg,uvw]=readbladed(filename)
% filename = 'C:\Users\ashim.giyanani\OneDrive - Windwise GmbH\Windwise\PostProcessing\TrialBladed_v2_writebl.wnd';

%eval(['fid=fopen(''',filename,'.wnd''',');']);
[Folder, baseFileName, extension] = fileparts(filename)
if isempty(extension)
  str = sprintf('%s.wnd',filename);
else
  str = filename
end  
fid = fopen(str);
New=fread(fid,1,'int16');
Spec=fread(fid,1,'uint16');

if Spec==3 | Spec==5 | Spec == 7
  Nturb=3;
else
  Nturb=1;
end
    
if Spec==4
  Nturb=fread(fid,1,'uint32');
  Latitude=fread(fid,1,'float32');
  z0=fread(fid,1,'float32');
  H=fread(fid,1,'float32');
  TI_u=fread(fid,1,'float32');
  TI_v=fread(fid,1,'float32');
  TI_w=fread(fid,1,'float32');
end

if Spec == 7
  unknown_1 = fread(fid, 1 , 'int16');
  unknown_2 = fread(fid, 1 , 'int16');
end

delta_z=fread(fid,1,'float32');
while delta_z<1
   delta_z = fread(fid, 1, 'float32');
end
delta_y=fread(fid,1,'float32');
delta_x=fread(fid,1,'float32');

Nx2=fread(fid,1,'uint32');
Nx=2*Nx2;

U=fread(fid,1,'float32');
Luz=fread(fid,1,'float32');
Luy=fread(fid,1,'float32');
Lux=fread(fid,1,'float32');

Dummy=fread(fid,1,'uint32');
Seed=fread(fid,1,'uint32');

Nz=fread(fid,1,'uint32');
Ny=fread(fid,1,'uint32');

x=[0:Nx-1]*delta_x;
% BLADED: first element at right bottom corner, i.e. maximum y
y=-[-(Ny-1)/2:(Ny-1)/2]*delta_y; % clockwise=false in Turbsim
z=[-(Nz-1)/2:(Nz-1)/2]*delta_z;
t=x/U;


if Nturb==3 | Spec==7
  Lvz=fread(fid,1,'float32');
  Lvy=fread(fid,1,'float32');
  Lvx=fread(fid,1,'float32'); % xLv
  Lwz=fread(fid,1,'float32');
  Lwy=fread(fid,1,'float32');
  Lwx=fread(fid,1,'float32'); % xLw
end

if Spec==7
   Alpha = fread(fid, 1, 'float32'); % coherence decay parameter
   Gamma = fread(fid, 1, 'float32'); % coherence scale parameter
end

u=nan(Nx,Ny,Nz);
v=nan(Nx,Ny,Nz);
w=nan(Nx,Ny,Nz);
Nplane=Ny*Nz*Nturb;
Ind1=1:3:Nplane;
Ind2=2:3:Nplane;
Ind3=3:3:Nplane;

for i=1:Nx,
   uvw=fread(fid,Nplane,'int16');
   if Nturb==3
     uu=uvw(Ind1);
     vv=uvw(Ind2);
     ww=uvw(Ind3);
     v(i,:,:)=reshape(vv,1,Ny,Nz);
     w(i,:,:)=reshape(ww,1,Ny,Nz);
   elseif Nturb==1
     uu=uvw;
   else
     error('wrong number of Nturb')  
   end
   u(i,:,:)=reshape(uu,1,Ny,Nz);
end

u=u./1000;
v=v./1000;
w=w./1000;

Clockwise = 0;
if Clockwise == 1
% rotate the u,v,w components to match turbsim output, look into Turbsim bladed format
u = rot90(permute(u,[3 2 1]),2);  % permute
v = rot90(permute(v,[3 2 1]),2);
w = rot90(permute(w,[3 2 1]),2);
y = flip(y); % min to max in left to right format
end;

if Spec == 4
    ug = U.*(TI_u.*0.01.*u+1); % Turbsim style U wind component
    vg = U.*(TI_v.*0.01.*v); % turbsim style V wind component
    wg = U.*(TI_w.*0.01.*w); % turbsim style W wind component
end
fclose(fid);
