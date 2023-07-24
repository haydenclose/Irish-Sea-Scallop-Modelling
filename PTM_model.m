%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% The model below is the particle tracking model developed by %%%
%%% Hayden Close, if you plan to use please contact me for any  %%%
%%% help at hclose_5@msn.com.                                   %%%
%%%                                                             %%%
%%% The model needs the oceanographic model developed by Peter  %%%
%%% Robins, p.robins@bangor.ac.uk, to run. In addition the      %%%
%%% generation of the random numbers (see folder) needs to be   %%%
%%% run before hand.                                            %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Clears everything
clear;                                                      % Clears all variables
clc;                                                        % Clears command window
clf;                                                        % Clears figure
fclose('all');                                              % Closes all files


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Generation of Random Numbers

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
np1= 60000;           % Iintial starting no of particles
tp = 10000;           % Number of particles used within the model
tt = 10272;           % Number of iterations (time steps)withing the model

rng(1);               % Uses a new seed, starting point to calculate random numbers
RN1=rand(np1,1);      % Postions for assigning the starting positon on the X direction
rng(2);
RN2=rand(np1,1);      % Postions for assigning the starting positon on the Y direction
rng(3);
RN3=rand(tp,1);       % For assigning the starting egg size
rng(4);
RN4=rand(tt,tp);      % For addition of a random element in the growth rate
rng(5);
RN5=rand(tt,tp);      % For use in determing the random walk (sub grid-scale turbulence)in the X direction
rng(6);
RN6=rand(tt,tp);      % For use in determing the random walk (sub grid-scale turbulence)in the Y direction
rng(7);
RN7=rand(tt,tp);      % For addition of a random element in the swiming speed at 0.5-1.0 mm s
rng(8);
RN8=rand(tt,tp);      % For addition of a random element in the swiming speed at 1.0-1.4 mm s
rng(9);
RN9=rand(tt,tp);      % For random element in the angle of the random walk



%% Load text files
celtic_sea = load('CelticSea.txt');                         % Load Celtic Sea coordinates
C_Bay = load('Cardigan_Bay.txt');                           % Load Cardigan Bay coordinates
Cornwall = load('North_Cornwall.txt');                      % Load North Cornwall coordinates
Tuskar = load('Tuskar.txt');                                % Load Tuskar coordinates
Trem_Bay = load('Tremadog Bay.txt');                        % Load Tremadog Bay coordinates
N_Wales = load('NorthWales.txt');                           % Load North Wales coordinates
IOM_Main = load('IOM mergercomplete.txt');                  % Load IoM main coordinates
IOM_Small =  load('top IoM patch.txt');                     % Load Point Ayre coordinates
Llyn = load ('Llyn Peninsula.txt');                         % Load Llyn Peninsula coordinates
SunriseSunset = load('SunRiseSunSet.txt');                  % Load Sunrise Sunset times, generated from the Suncycle.m script (available online)
% load('RandNos.mat');                                        % Loads the random numbers needed for the model (Needs to be run code first)
% key = 1:10000;                                              % Creates a key for the calculation of the mortality

%% Find and load the grid and bathymetry, only need to generate once                   
% folder = 'U:\College of Natural Sciences\Research Data\CAMS-GROUP-PeterRobins\mar1990\'; % Data on U drive
folder = 'E:\PTM model data\';                              % My hardrive day data
depthandgrid = ([folder,'mar90.00001.nc']);                 % Selects first file
ncid = netcdf.open(depthandgrid,'NC_NOWRITE');              % Opens file
depth = netcdf.getVar(ncid,0,'double');                     % Loads the bathymetry
landsea = netcdf.getVar(ncid,1,'double');                   % Loads the border between land and sea
netcdf.close(ncid);                                         % Closes netcdf to make MATLAB quicker
landsea(landsea==0)=NaN;                                    % replaces 0 with NaN
depth(depth==0)=NaN;                                        % replaces 0 with NaN

%% Irish Sea and POM variables
X1 = -9.0; X2 = -2.7;                                       % X-axis
Y1 = 50.1; Y2 = 56.;                                        % Y-axis
xspac = 1/30;                                               % X-axis grid distances
yspac = 1/60;                                               % Y-axis grid distances
long = X1: xspac: X2;                                       % Long = x axis
lati = Y1: yspac: 55.99;                                    % Lat = y axis 55.99 (slight adjustment due to curvature of Earth)
nt   = 10272;                                               % Number of time steps March-September
nx = 190; ny = 354; nz = 21;                                % Grid dimensions of the model
dt = 1800;                                                  % POM output time step, 1800 seconds (30 mins)
dx = 2000.; dy = 2000.;                                     % Cell size (m), (2km by 2km)
year = 1990;                                                % Year of model
%% Records movie
record = 0;                                                 % Record model (1) or not
if record ==1;                                              % Record model (1) or not
vid = VideoWriter('April Cardigan Bay.avi');                % Creates Movie, change name of movie here
vid.Quality = 100;                                          % Quality of video 1-100 (100 = best)
vid.FrameRate = 146;                                        % No frames per Second
open(vid);                                                  % Opens the video for recording
end                                                         % Closes record movie if
%%
polygon = 1;                                                % 1 = plot particles in desired polygon defined on line 83 & 84, 0 = plot particles outside polygons
avoidland = 1;                                              % 1 = plot particles avoiding land even if polygon 0 = plots particles on land
start = 1488;                                               % Starting iteration defines date
finish = 2000;                                              % Finishing iteration
np1 = 50000;                                                  % Number of particles released
idim = 2;                                                   % 2= plot in 2D , 3= plot in 3D variables

%% Shape files Pre allocation for speed
xpCeltic=zeros(216,1);   ypCeltic=zeros(216,1);             % Celtic Patch
xpC_Bay=zeros(34,1);     ypC_Bay=zeros(34,1);               % Cardigan Bay Patch
xpTuskar=zeros(53,1);    ypTuskar=zeros(53,1);              % Tuskar Patch
xpIOM_M=zeros(537,1);    ypIOM_M=zeros(537,1);              % Main IoM Patch
xpIOM_S=zeros(36,1);    ypIOM_S=zeros(36,1);                % Point Ayre Patch
xpNW=zeros(20,1);        ypNW=zeros(20,1);                  % North Wales Patch
xpTB=zeros(9,1);         ypTB=zeros(9,1);                   % Tremadog Bay Patch
xpCwall=zeros(51,1);     ypCwall=zeros(51,1);               % North Corwall Patch
xpLlyn=zeros(8,1);       ypLlyn=zeros(8,1);                 % Llyn Peninsula Patch

%% Defines the X and Y coordiinates for each patch
for CS = 1:216;    xpCeltic(CS)= celtic_sea(CS,1);      ypCeltic(CS)= celtic_sea(CS,end);     end; % Celtic Patch                                       
for CBay = 1:34;   xpC_Bay(CBay)= C_Bay(CBay,1);        ypC_Bay(CBay)= C_Bay(CBay,end);       end; % Cardigan Bay Patch                                                        
for Tusk = 1:53;   xpTuskar(Tusk)= Tuskar(Tusk,1);      ypTuskar(Tusk)= Tuskar(Tusk,end);     end; % Tuskar Patch 
for IOM_M = 1:537; xpIOM_M(IOM_M)= IOM_Main(IOM_M,1);  ypIOM_M(IOM_M)= IOM_Main(IOM_M,end);   end; % Main IoM Patch  
for IOM_S = 1:36;  xpIOM_S(IOM_S)= IOM_Small(IOM_S,1);  ypIOM_S(IOM_S)= IOM_Small(IOM_S,end); end; % Point Ayre Patch
for NW = 1:20;     xpNW(NW)= N_Wales(NW,1);             ypNW(NW)= N_Wales(NW,end);            end; % North Wales Patch
for TB = 1:9;      xpTB(TB)= Trem_Bay(TB,1);            ypTB(TB)= Trem_Bay(TB,end);           end; % Tremadog Bay Patch
for CW = 1:51;     xpCwall(CW)= Cornwall(CW,1);         ypCwall(CW)= Cornwall(CW,end);        end; % North Corwall Patch
for Ll = 1:8;      xpLlyn(Ll)= Llyn(Ll,1);               ypLlyn(Ll)= Llyn(Ll,end);            end; % Llyn Peninsula Patch

%% Initial example release site:
%IOMish
xmin = 109*dx; xmax = 143*dx; ymin = 117*dy; ymax = 155*dy; zmin= 21; zmax = 21; % Gives  start points in metres
polyXS = xpC_Bay;                                                                % Defines start polygon for X
polyYS = ypC_Bay;                                                                % Defines start polygon for Y
xFi = zeros(np1,1); yFi = zeros(np1,1);                                          % Pre-defines matrix with zeros for x, y and z to make matlab quicker

%% Loops for all particles to give a random position in a polygon
npcount=0;                                                    % Start for count to 10,000 so the number of particles always 10,000
for Fi=1:np1;                                                 % Loop for all First Intial particles
xFi(Fi)=(RN1(Fi)*(xmax-xmin))+xmin;                           % Random place in the grid with the x direction
yFi(Fi)=(RN2(Fi)*(ymax-ymin))+ymin;                           % Random place in the grid with the y direction
xx=X1+(xFi/dx)*xspac;                                         % converts x into coordinates to compare against polygons
yy=Y1+(yFi/dy)*yspac;                                         % Converts y into coordinates to compare against polygons
    if polygon ==1;                                           % If polygon =1 then particles limited =1 to the defined patch
        TB_Boundary = inpolygon(xx,yy,polyXS,polyYS);         % Defines whether points in the polygon 1 = in boundary 0 = outside
        if   TB_Boundary(Fi) == 1;                            % If 1 (in polygon) keeps x value, if outside polygon x given a value of 0
        xFi(Fi) = xFi(Fi); yFi(Fi) = yFi(Fi);                 % x(i) = x, y = y 
        else xFi(Fi)=0; yFi(Fi)=0;                            % Else x(i) = 0 and y = 0
        end;                                                  % End, give x a value of 0 if outside of polygon
    else xFi(Fi) = xFi(Fi); yFi(Fi) = yFi(Fi);                % if polygon doesnt = 1 all x keep x value regardless of wether inside our outside of polygon
    end ;                                                     % End if inside or outside of particle
end;                                                          % End random place for particle
xin=xFi(xFi~=0);                                              % Removes the 0s for x given to those outside of polygon
yin=yFi(yFi~=0);                                              % Removes the 0s for y given to those outside of polygon
    if   polygon == 1;                                        % If plotting within polygon                                                             
         LS=sum(TB_Boundary);                                 % New no of particles = total inside of defined polygon
    else LS=Fi;                                               % If polygon does not = 1 then LS = same as Fi
    end;                                                      % End if for creation of new no of particles
x = zeros(LS,1);                                              % Preload matrix to make MATLAB quicker
y = zeros(LS,1);                                              % Preload matrix to make MATLAB quicker

 %% stop particles being located on land
for ii=1:LS                                                   % Loop from 1 to new no of particles to define I and J to determine if particle given a postion on land
II = floor(xin(ii)/dx)+1;                                     % Locates cell in x direction in m, +1 is to take account of the lack of an 0 cell
JJ = floor(yin(ii)/dy)+1;                                     % Locates cell in y direction in m, +1 is to take account of the lack of an 0 cell
     if avoidland == 1;                                       % If avoidland =1 then particles avoid land
         if   isnan(landsea(II,JJ))==0;                       % Isnan = 0 then particles in sea if isnan = 1 then particles on land
              x(ii) = xin(ii); y(ii) = yin(ii);               % x = x, y = y 
         else x(ii)=0; y(ii)=0;                               % If particles randomly placed on land x = 0, y = 0
        end;                                                  % End if isnan = 0 for on land
    else x(ii) = xin(ii); y(ii) = yin(ii);                    % If dont want to avoid land all x's= x all y's = y 
    end;                                                      % End avoidland = 1
end;                                                          % end locate postion 
x=x(x~=0);                                                    % Removes the 0s for x
y=y(y~=0);                                                    % Removes the 0s for y
np=numel(x);                                                  % New no of particles
for lim = 1:np;                                               % Loops through new no of particles
if npcount <10000                                             % Limits to 10000 only      
    npcount = npcount +1;                                     % Counts to 10000
    x(lim) = x(lim); y(lim) = y(lim);                         % If under 10000 then x = x 
else x(lim)= 0; y(lim)= 0;                                    % If over 10000 then x = 0
end                                                           % End loop
end                                                           % End loop                          
x=x(x~=0);                                                    % Removes the 0s for x
y=y(y~=0);                                                    % Removes the os for y
np=numel(x);                                                  % Calculates number of particles

%% Finding position and depth and starting size of larvae
z =  zeros(np,1);                                             % Predefining Z to make MATLAB quicker
larvaeLength = zeros(np,1);                                   % Predefining larvaeLength to make MATLAB quicker            
for iii=1:np;                                                 % Loop to determine position of particles in I and J of particles                                      
    I = floor(x(iii)/dx)+1;                                   % Locates cell in x direction in m, +1 is to take account of the lack of an 0 cell
    J = floor(y(iii)/dy)+1;                                   % Locates cell in y direction in m, +1 is to take account of the lack of an 0 cell
    K = 20;                                                   % Starting depth (sea floor)
    z(iii)=depth(I,J);                                        % Depth at the location I,J
% Ititial starting size of egg
    eggmin = 64.2;                                            % Min egg size
    eggmax = 72.2;                                            % Max egg size
   larvaeLength(iii) = eggmin + (RN3(iii)*(eggmax-eggmin));   % Gives each particle a random egg size between to min and max
end;                                                          % End locate new position and start larvae length

%% Predefining matrixes to make MATLAB quicker
settle=zeros(np,1);                                           % % Pre-defines Settle matrix with zeros
stage = ones(np,1);                                           % Predefines particle status (0=dead; 1-5=alive, stages = 1-5), larvae start as 1 as eggs
u = zeros(np,1); v = zeros(np,1); w = zeros(np,1);            % Pre-defines U V W matrixes with zeros
T = zeros(np,1);                                              % Pre-defines Temperature matrix with zeros
swim = zeros(np,1);                                           % Pre-defines swiming as zero, early stages cannot swim
growthrate = zeros (np,1);                                    % Pre-defines growthrate matrix with zeros 
ddx_sq = zeros (np,1); ddy_sq = zeros (np,1); dist = zeros (np,1); % Pre-defines distance matrixes with zeros 
L_die= zeros (np,1);                                          % Pre-defines length at which larvae die matrix with zeros 
distance = zeros (np,1);                                      % Pre-defines distance matrix with zeros 
Rec_growthrate = zeros (finish,np);                           % Pre-defines record growth rate matrix with zeros 
Rec_zprofile = zeros (finish,np);                             % Pre-defines record depth profile matrix with zeros 
Rec_stage = zeros (finish,np);                                % Pre-defines record stage matrix with zeros 
Rec_t_distance = zeros (finish,np);                           % Pre-defines record total distance matrix with zeros
Rec_temperature = zeros (finish,np);                          % Pre-defines record temperature matrix with zeros 
Rec_settleday = zeros (finish,np);                            % Pre-defines record settlementday matrix with zeros 
Rec_x = zeros (finish,np);                                    % Pre-defines record x matrix with zeros 
Rec_y = zeros (finish,np);                                    % Pre-defines record y matrix with zeros 
%% MAIN LOOP inital code
Filelist=dir([folder,'mar90.*']);                             % Finds the files for 1:x with PTM data
Filelisttemp=dir([folder,'restart.*']);                       % Finds the files for 1:x for other variables (temperature only used)
% Load in initial temperature matricies: T0
itemp=floor(start/96)+1;                                      % Find the initial file number to read, (96 iterations every two days)
filetemp=([folder,Filelisttemp(itemp).name]);                 % Select initial restart file (contains environmental data)
ncid = netcdf.open(filetemp,'NC_NOWRITE');                    % Opens initial restart file(contains environmental data)
T0 = netcdf.getVar(ncid,24,'double');                         % Initial temperature to expolate between as only got data every 2 days
netcdf.close(ncid);                                           % Closes netcdf to quicken MATLAB up
% T1
filetemp=([folder,Filelisttemp(itemp+1).name]);               % Select initial restart file (contains environmental data)
ncid = netcdf.open(filetemp,'NC_NOWRITE');                    % Opens initial restart file (contains environmental data)
T1 = netcdf.getVar(ncid,24,'double');                         % Second temperature to extrapolate to
netcdf.close(ncid);                                           % Closes netcdf to quicken MATLAB up

for it=start:finish;
    %% Calculations for the month, day and hour of each iteration
  if it<1488;              month = 03; MonthName = 'March';     day=fix(it/48)+1;        hour = (it-((day-1)*48))/2; end       % Month 3 = March
  if it>=1488 && it<2928;  month = 04; MonthName = 'April';     day=fix((it-1488)/48)+1; hour = (it-1488-((day-1)*48))/2; end  % Month 4 = April
  if it>=2928 && it<4416;  month = 05; MonthName = 'May';       day=fix((it-2928)/48)+1; hour = (it-2928-((day-1)*48))/2; end  % Month 5 = May
  if it>=4416 && it<5856;  month = 06; MonthName = 'June';      day=fix((it-4416)/48)+1; hour = (it-4416-((day-1)*48))/2; end  % Month 6 = June
  if it>=5856 && it<7344;  month = 07; MonthName = 'July';      day=fix((it-5856)/48)+1; hour = (it-5856-((day-1)*48))/2; end  % Month 7 = July
  if it>=7344 && it<8832;  month = 08; MonthName = 'August';    day=fix((it-7344)/48)+1; hour = (it-7344-((day-1)*48))/2; end  % Month 8 = August
  if it>=8832 && it<10272; month = 09; MonthName = 'September'; day=fix((it-8832)/48)+1; hour = (it-8832-((day-1)*48))/2; end  % Month 9 = September

file=([folder,Filelist(it).name]);                            % Select hydrodynamic file
ncid = netcdf.open(file,'NC_NOWRITE');                        % Opens file
if idim==2;                                                   % If using 2D matrix to test model works
    U0 = netcdf.getVar(ncid,2,'double');                      % U velocity:(U averaged)
    V0 = netcdf.getVar(ncid,3,'double');                      % V velocity:(V averaged)
    [U]=u2rho_2d(U0);                                         % Interpolate to 'centered' grid (Need the function in the MATLAB folder)
    [V]=v2rho_2d(V0);                                         % Interpolate to 'centered' grid (Need the function in the MATLAB folder)
else                                                          % For when wanting to use 3D matrix
    U0 = netcdf.getVar(ncid,5,'double');                      % U velocity
    V0 = netcdf.getVar(ncid,6,'double');                      % V velocity
    W = netcdf.getVar(ncid,7,'double');                       % W (vertical) velocity: 3D matrix
    [U]=u2rho_2d(U0(:,:,K));                                  % Interpolate to 'centered' grid (Need the function in the MATLAB folder)
    [V]=v2rho_2d(V0(:,:,K));                                  % Interpolate to 'centered' grid (Need the function in the MATLAB folder)
end                                                           % End for 2D or 3D      
 
% update temperature data every 2 days:          
if mod(it,96)==0 && itemp<106;                                % every two day
    itemp=itemp+1;                                            % Next iteration to interpolate between the two values
    T0=T1;                                                    % T1 now becomes zero
    filetemp=([folder,Filelisttemp(itemp+1).name]);           % Select restart file
    ncid = netcdf.open(filetemp,'NC_NOWRITE');                % Opens restart file
    T1 = netcdf.getVar(ncid,24,'double');                     % Temperature: 3D matrix
      netcdf.close(ncid);                                     % Closes netcdf to quicken matlab
end                                                           % End update temperature
temp = T0 + ((T1-T0)/(2*24*3600))*(mod(it,96)*dt);            % Interpolates temperature to local time

%% Loop for np particles
for i=1:np;                                                   % Loop over particles:
    I = floor(x(i)/dx)+1;                                     % Locates cell in x direction in m, +1 is to take account of the 0
    J = floor(y(i)/dy)+1;                                     % Locates cell in y direction in m, +1 is to take account of the 0
    u(i) = U(I,J) + ((U(I+1,J)-U(I,J))/dx)*((x(i)-(I-1)*dx)); % Velocity U in I and J direction
    v(i) = V(I,J) + ((V(I,J+1)-V(I,J))/dx)*((y(i)-(J-1)*dy)); % Velocity V in I and J direction
if idim==2;                                                   % 2D 
    T(i) = temp(I,J,1);                                       % Temperature for position I,J and K(depth)
else                                                          % Run 3D model parameters for depth loacation
    K = floor( (z(i)/depth(I,J))*(nz-1) );                    % Locates depth
    if(K<1);     K=1;    end                                  % Cant go out of the sea vertically into the air
    if(K>nz-1);  K=nz-1; end                                  % Cant go out of the sea vertically into the seabed
    w(i) = W(I,J,K);                                          % Vertical velocity
    T(i) = temp(I,J,K);                                       % Temperature for position I,J and K(depth)
end                                                           % End 3D model for depth location
if T(i)> 17; T(i)= 17; end;                                   % Temp cant go over 17 C

 %% Growth Rate
if T(i)  <= 9.5;                                              % Temp range 8.5-9.5
    lowlimit = 1.36;                                          % Lower growth rate
    range  = 0.5;                                             % Difference between upper and lower values
    growthrate(i) = (lowlimit + RN4(it,i)*range)/(24*3600/dt);% Calculates growthrate
    elseif T(i) > 9.5 && T(i) <= 10.5;                        % Temp range 9.5-10.5
    lowlimit = 1.86;                                          % Lower growth rate
    range  = 0.5;                                             % Difference between upper and lower values
    growthrate(i) = (lowlimit + RN4(it,i)*range)/(24*3600/dt);% Calculates growthrate
    elseif T(i) > 10.5 && T(i) <= 11.5;                       % Temp range 10.5-11.5
    lowlimit = 2.36;                                          % Lower growth rate
    range  = 0.5;                                             % Difference between upper and lower values
    growthrate(i) = (lowlimit + RN4(it,i)*range)/(24*3600/dt);% Calculates growthrate
    elseif T(i) > 11.5 && T(i) <= 12.5;                       % Temp range 11.5-12.5
    lowlimit = 2.87;                                          % Lower growth rate
    range  = 0.52;                                            % Difference between upper and lower values
    growthrate(i) = (lowlimit + RN4(it,i)*range)/(24*3600/dt);% Calculates growthrate
    elseif T(i) > 12.5 && T(i) <= 13.5;                       % Temp range 12.5-13.5
    lowlimit = 3.29;                                          % Lower growth rate
    range  = 0.31;                                            % Difference between upper and lower values
    growthrate(i) = (lowlimit + RN4(it,i)*range)/(24*3600/dt);% Calculates growthrate
    elseif T(i) > 13.5 && T(i) <= 14.5;                       % Temp range 13.5-14.5
    lowlimit = 3.60;                                          % Lower growth rate
    range  = 0.31;                                            % Difference between upper and lower values% Lower growth rate
    growthrate(i) = (lowlimit + RN4(it,i)*range)/(24*3600/dt);% Calculates growthrate
    elseif T(i) > 14.5 && T(i) <= 15.5;                       % Temp range 14.5-15.5
    lowlimit = 3.91;                                          % Lower growth rate
    range  = 0.32;                                            % Difference between upper and lower values
    growthrate(i) = (lowlimit + RN4(it,i)*range)/(24*3600/dt);% Calculates growthrate
    elseif T(i) > 15.5 && T(i) <= 16.5;                       % Temp range 15.5-16.5
    lowlimit = 4.26;                                          % Lower growth rate
    range  = 0.38;                                            % Difference between upper and lower values
    growthrate(i) = (lowlimit + RN4(it,i)*range)/(24*3600/dt);% Calculates growthrate
    elseif T(i) > 16.5 && T(i) <= 17.5;                       % Temp range 16.5.5-17.5
    lowlimit = 4.64;                                          % Lower growth rate
    range  = 0.38;                                            % Difference between upper and lower values
    growthrate(i) = (lowlimit + RN4(it,i)*range)/(24*3600/dt);% Calculates growthrate
    elseif T(i) > 17.5 && T(i) <= 18.5;                       % Temp range 17.5-18.5
    lowlimit = 5.02;                                          % Lower growth rate
    range  = 0.38;                                            % Difference between upper and lower values
    growthrate(i) = (lowlimit + RN4(it,i)*range)/(24*3600/dt);% Calculates growthrate
end                                                           % End temperature / growth rate relationship
% Calulation of new length of the larvae
if settle (i) == 0;                                           % Calculates new larval length if not yet settled
    larvaeLength(i) = larvaeLength(i)+ growthrate(i);         % New larval length
else larvaeLength(i) = larvaeLength(i);                       % If settled keeps same larval length
end                                                           % ends calculate new length
% Larval stage designation    
if stage(i) == 0;                                             % If larvae dead, stage = 0
   stage(i) =0;                                               % larvae keep dead status
elseif stage(i) == 7;                                         % if larvae left domain, stage = 7
   stage(i) =7;                                               % Larve stay seperate
     elseif larvaeLength(i) >= 64.2 && larvaeLength(i) < 76;  % Larvae length range 64.2 - 76
     stage(i) = 1;                                            % Egg to Trochophore
     elseif larvaeLength(i) >= 77 && larvaeLength(i) < 81;    % Larvae length range 77 - 81
     stage(i) = 2;                                            % Trocophore and early veliger
     elseif larvaeLength(i) >= 82 && larvaeLength(i) < 183;   % Larvae length range 82 - 183
     stage(i) = 3;                                            % Veliger
     elseif larvaeLength(i) >= 184 && larvaeLength(i) < 214;  % Larvae length range 184 - 214
     stage(i) = 4;                                            % Eyed veliger  
     elseif larvaeLength(i) >= 215 && larvaeLength(i) < 240;  % Larvae length range 215 - 240
    stage(i) = 5;                                             % Pediveliger
           elseif larvaeLength(i) >= 241;                     % After 241 larvae settle
     stage(i) = 6;                                            % If pediveligers and not settled 
     L_die(i) = L_die(i) + 1;                                 % Counts no of days since able to spawn
end                                                           % End stage designition
if L_die(i)  >= 720 && stage(i) ~=7;                          % If no of days greater than 14 and stage not = to 7
stage(i) = 0;                                                 % Larvae 'die'
end                                                           % End delay period if

 %% Movement of Particles
% Horizontal Diffusion (random walk)
r = 1/sqrt(6);                                                % Standard deviation
KD  = hypot(U(i),V(i))/100;                                   % Horizontal diffusion from POM output (hypot robust computation of sqaure root of the sum of squares)
ang = 2*pi*RN9(it,i);                                         % Angle of the random movement
Lx = (RN5(it,i)/r).*(cos(ang)).*(sqrt(2*KD*dt));              % Movement in the x direction
Ly = (RN6(it,i)/r).*(sin(ang)).*(sqrt(2*KD*dt));              % Movement in the y direction
ddx = (u(i)*dt)+Lx; ddy = (v(i)*dt)+Ly;                       % New movement in x and y directions
x(i) = x(i)+ddx;                                              % Position from previous iteration plus new movement in x direction
y(i) = y(i)+ddy;                                              % Position from previous iteration plus new movement in x direction 
dist(i) = sqrt(ddx^2+ddy^2);                                  % distance = sqrt(x^2 + y^2)
distance(i) = distance(i) + dist(i);                          % Total distance = previous distance + distance travelled in interation

%% 3D to calculate swimming speed and Diel migrations
if idim==3;                                                   % 3D model
% swimming speed:
    if     stage(i) == 0;                                     % If stage = 0 then larvae dead
           swimspd = 0;                                       % cannot swim
    elseif stage(i) == 1;                                     % If stage = 0 then larvae dead
           swimspd = 0;                                       % cannot swim
    elseif stage(i) > 2;                                      % Trochophore to early veliger
           swimspd = ((0.5 + 0.5*RN7(it,i))/1000)*1800;       % Swim in range 0.5 to 1   
    else   swimspd = ((1 + 0.4*RN8(it,i))/1000)*1800;         % Swim in range 1 to 1.4
    end;                                                      % End if for if 3D so swimming speed can be included

%% Vertical diel Migrations
% call the_suncycle function:
date = [year month day];                                      % Calculates date
    if hour>SunriseSunset(it,1) & hour<SunriseSunset(it,end) & stage(i)==3;    % If hour within day light hours and the individual larva a veliger
    swim(i) = -swimspd;                                       % Day: Swim = -ve
    else  swim(i) = +swimspd;                                 % If night, larvae swim up
    end;                                                      % End determine swim direction 
if stage (i)> 3                                               % If eyed-veligers or later 
    swim(i) = -swim(i);                                       % Swim to sea bed
end                                                           % End swim if
ddz = w(i)*dt;                                                % Vertical movement
z(i) = z(i)-ddz-(swim(i));                                    % Movement vertically
if stage(i) == 3 & z(i) > 10;                                 % If veliger and max depth exceeds 10 m
   z(i) = 10;                                                 % Depth = 10 m
else                                                          % else
end;                                                          % End limit depth
end;                                                          % End if 3d model to calculate swim

%% Check particles make sure dont leave the domain and kills them
% Check no particles leave horizontal domain:
if x(i)<= 1;     x(i)=1; stage(i)=7;  end;                    % If less than the 1st grid (cant be in 0 grid) puts back to 1 and killed
if x(i)>= nx*dx; x(i)=1; stage(i)=7;  end;                    % If greater than the 1st grid (cant be in 0 grid) taken to 1 and killed
if(y(i)<= 1);    y(i)=1; stage(i)=7;  end;                    % If less than the 1st grid (cant be in 0 grid) puts back to 1 and killed
if(y(i)>= ny*dy);y(i)=1; stage(i)=7;  end;                    % If greater than the 1st grid (cant be in 0 grid) taken to 1 and killed
% Check no particles leave vertical domain:
if(z(i)<0);          z(i)=0;          end;                    % If x < 0 (0m = top) then depth = 0 keeps at surface water
if(z(i)>depth(I,J)); z(i)=depth(I,J); end;                    % If z > than depth of sea then z = maximum depth

%% Reflection of particles
I = floor(x(i)/dx)+1;                                         % Locates cell in x direction in m, +1 is to take account of the 0
J = floor(y(i)/dy)+1;                                         % Locates cell in y direction in m, +1 is to take account of the 0
land = isnan(landsea(I,J));                                   % If particle hit land 1 = land, NaN = sea
if land == 1;                                                 % Then particle has hit land
     x(i)= x(i)-ddx;                                          % Reflect back to previous position
     y(i)= y(i)-ddy;                                          % Reflect back to previous position
end;                                                          % End Relection of particles loop
end                                                           % Ends particle loop (np)
%% Record variables
xUzero=X1+(x/dx)*xspac; yUzero=Y1+(y/dy)*yspac;               % Converts X an Y into coordinates
Rec_x(it,:)= xUzero;                                          % Records X for every larvae for every iteration
Rec_y(it,:)= yUzero;                                          % Records Y for every larvae for every iteration                
Rec_growthrate(it,:)= larvaeLength;                           % Records growthrate for every larvae for every iteration
Rec_stage(it,:) = stage;                                      % Records stage for every larvae for every iteration
Rec_t_distance (it,:) = distance;                             % Records total distance for every larvae for every iteration
Rec_temperature (it,:) = T;                                   % Records record temperature for every larvae for every iteration
Rec_zprofile(it,:)= z;                                        % Records depth profile for every larvae for every iteration
minutes = 30*(it-start);                                      % Converts the iteration into minutes
hour = minutes/60;                                            % Converts minutes into hours
Day = (hour/24);                                              % Converts hours to days
dayoutput= fix(Day);                                          % Supresses the output of the decimal places

%% Plot:
clf;                                                          % Clear figure so dont get overlapping figures
hold on;                                                      % Hold on allows multiple additions to same figure without over writing
sst=temp(:,:,1);                                              % Get temp for I,J,K
for ii=1:nx;                                                  % Loop II
    for jj=1:ny;                                              % Loop JJ
        if sst(ii,jj) > 17; sst(ii,jj) = 17; end;             % If sst > 17 then it is limited to 17
    end;                                                      % End jj
end;                                                          % End ii
PTM = pcolor(long, lati, landsea'.*sst');                     % Sea seasurface temperature and plot of the land
set(PTM, 'EdgeColor', 'none');                                % Removes boxes
shading flat                                                  % Shading of Plot
xp=X1+(x/dx)*xspac;                                           % Converts X into a coordinate position
yp=Y1+(y/dy)*yspac;                                           % Converts Y into a coordinate position
Celticpatch=plot(xpCeltic,ypCeltic,'k');                      % Plots Celtic scallop patch
CBaypatch=plot(xpC_Bay,ypC_Bay,'k');                          % Plots Cardigan Bay scallop patch
Tuskarpatch=plot(xpTuskar,ypTuskar,'k');                      % Plots Tuskar scallop patch
Trem_Baypatch=plot(xpTB,ypTB,'k');                            % Plots Tremadog Bay scallop patch
IOM_Mpatch=plot(xpIOM_M,ypIOM_M,'k');                         % Plots IoM main scallop patch
IOM_Spatch=plot(xpIOM_S,ypIOM_S,'k');                         % Plots IoM top small scallop patch
NWpatch=plot(xpNW,ypNW,'k');                                  % Plots North Wales scallop patch
Cornwallpatch=plot(xpCwall,ypCwall,'k');                      % Plots Cornwall Patch scallop patch
Llynpatch=plot(xpLlyn,ypLlyn,'k');                            % Plots Llyn peninsula scallop patch
a=find(stage(:)==0); plot((X1+(x(a)/dx)*xspac), (Y1+(y(a)/dy)*yspac),'.w');     % Plots stage 0(dead)particles 
b=find(stage(:)==1); plot((X1+(x(b)/dx)*xspac), (Y1+(y(b)/dy)*yspac),'.k');     % Plots stage 1(egg-grastula) particles
c=find(stage(:)==2); plot((X1+(x(c)/dx)*xspac), (Y1+(y(c)/dy)*yspac),'.y');     % Plots stage 2(trocophore-early veliger) particles
d=find(stage(:)==3); plot((X1+(x(d)/dx)*xspac), (Y1+(y(d)/dy)*yspac),'.b');     % Plots stage 3(Veligers)particles
e=find(stage(:)==4); plot((X1+(x(e)/dx)*xspac), (Y1+(y(e)/dy)*yspac),'.g');     % Plots stage 4(eyed-veligers)particles
f=find(stage(:)==5); plot((X1+(x(f)/dx)*xspac), (Y1+(y(f)/dy)*yspac),'.m');     % Plots stage 5(pediveligers)particles
g=find(stage(:)==6); plot((X1+(x(g)/dx)*xspac), (Y1+(y(g)/dy)*yspac),'.r');     % Plots stage 6(able to settl)particles
h=find(stage(:)==7); plot((X1+(x(h)/dx)*xspac), (Y1+(y(h)/dy)*yspac),'.k');     % Plots stage 7(left domain)particles
set(gca,'color', [0 .6 0],...                                 % Changes background colour(map) to green (land)
    'DataAspectRatio',[1 yspac/xspac 1]);                     % Plots axis properly so not stretched 
%axis([-5.75 -4 52 53]);                                      % Plot of a specific area in figure
axis([X1 X2 Y1 Y2]);                                          % Axis of full extent of the Irish Sea region
cbar = colorbar('location', 'southoutside');                  % Location  of colour bar
caxis([6 17]);                                                % Colour bar axis
set(cbar, 'YColor', [0 0 0]);                                 % Set colour bar text to black, x
set(cbar, 'XColor', [0 0 0]);                                 % Set colour bar text to black, y
xlabel ('Longitude');                                         % x-axis label
ylabel('Latitude');                                           % y-axis label
title( 'Pecten maximus Larval Dispersion Model');             % Title of plot
text(-9.25, 56.5, ['Larval Duration = ' num2str(dayoutput)]); % Position and display of the day on the figure
text(-5.75, 56.5, ['Date = ' num2str(day),' ',num2str(MonthName),' ',num2str(year)],...    
    'VerticalAlignment','middle',...
    'HorizontalAlignment','center');                          % Plots date on figure
text(-2, 56.5, [' Iteration = ' num2str(it)]);                % Plots iteration number on figure
if   record == 1;                                             % If recording movie
     writeVideo(vid, getframe(gcf));                          % Get frame assembles array of frames to playback using WriteVideo
else getframe(gcf);                                           % Get frame
end                                                           % End get frame
end                                                           % End (it) loop

%% calculate position of larvae when settled
% finds when the larvae were in the sites
Celtic_Inpolgon = inpolygon(Rec_x,Rec_y,xpCeltic,ypCeltic);   % Output for Celtic scallop patch
C_Bay_Inpolygon = inpolygon(Rec_x,Rec_y,xpC_Bay,ypC_Bay);     % Output for Cardigan Bay scallop patch
Tuskar_Inpolygon = inpolygon(Rec_x,Rec_y,xpTuskar,ypTuskar);  % Output for Tuskar scallop patch
Trem_Bay_Inpolygon = inpolygon(Rec_x,Rec_y,xpTB,ypTB);        % Output for Tremadog Bay scallop patch
IOM_main_Inpolygon = inpolygon(Rec_x,Rec_y,xpIOM_M,ypIOM_M);  % Output for IoM main scallop patch
IOM_small_Inpolygon = inpolygon(Rec_x,Rec_y,xpIOM_S,ypIOM_S); % Output for IoM top small scallop patch
N_Wales_Inpolygon = inpolygon(Rec_x, Rec_y,xpNW,ypNW);        % Output for North Wales scallop patch
Cornwall_Inpolygon = inpolygon(Rec_x,Rec_y,xpCwall,ypCwall);  % Output for Cornwall Patch scallop patch
Llyn_Inpolygon = inpolygon(Rec_x,Rec_y,xpLlyn,ypLlyn);        % Output for Llyn peninsula scallop patch
settle = zeros (finish,np);                                   % Pre-defining settle to make MATLAB quicker
OUTx = zeros (np,1);                                          % Pre-defining Xoutput to make MATLAB quicker
OUTy = zeros (np,1);                                          % Pre-defining Youtput to make MATLAB quicker
OUTiteration = zeros (np,1);                                  % Pre-defining interation output to make MATLAB quicker
OUTxa = zeros (np,1);                                         % Pre-defining alive larve's X coordinate to make MATLAB quicker
OUTya = zeros (np,1);                                         % Pre-defining alive larve's Y coordinate to make MATLAB quicker
OUTDistance  = zeros (np,1);                                  % Pre-defining distance output to make MATLAB quicker
OUTgrowth = zeros (np,1);                                     % Pre-defining growth output to make MATLAB quicker
OUTxd = zeros (np,1);                                         % Pre-defining dead larve's X coordinate to make MATLAB quicker
OUTyd = zeros (np,1);                                         % Pre-defining dead larve's Y coordinate to make MATLAB quicker
OUTsettled = zeros (np,1);                                    % Pre-defining settled output to make MATLAB quicker
OUTstage = zeros (np,1);                                      % Pre-defining stage output to make MATLAB quicker
for it = start:finish;                                        % Loop from start to finish
for i = 1:10000                                               % Loop for all particles 
if C_Bay_Inpolygon(it,i) ==1 || IOM_main_Inpolygon(it,i) ==1 || Celtic_Inpolgon(it,i) ==1 || Tuskar_Inpolygon(it,i) ==1 || Trem_Bay_Inpolygon(it,i)==1 || IOM_main_Inpolygon(it,i) ==1 || IOM_small_Inpolygon(it,i) ==1 || N_Wales_Inpolygon(it,i) ==1 || Cornwall_Inpolygon(it,i) ==1  || Llyn_Inpolygon (it,i)==1;
settle(it,i)=1;                                               % If particle inside settle = 1
end                                                           % end if inpolygon
end                                                           % end loop for it
end                                                           % end loop for i 
Rec_stage2 = Rec_stage;                                       % Create copy of record stage
for it = 1488:finish;                                         % Loop for it 
for i = 1:10000;                                              % Loop for i
if Rec_stage2 (it, i)== 7;                                    % If stage = 7, out of domain
   Rec_stage2(it, i) = 0;                                     % stage = 0 'dead' larvae
end                                                           % end if
if Rec_stage2 (it,i)== 6 && settle(it,i)== 1 && OUTgrowth(i)==0;% If stage = 6 (able to settle), settle = 1 (settled) and not already been assigned
OUTx(i) = Rec_x(it,i);                                        % Records the X coordinate for time of settlement for all the larvae
OUTy(i) = Rec_y(it,i);                                        % Records the Y coordinate for time of settlement for all the larvae
OUTiteration (i)= it;                                         % Records the iteration output for time of settlement for all the larvae
OUTxa(i) = Rec_x(it,i);                                       % Records the X coordinate for time of settlement for all the alive larvae
OUTya(i) = Rec_y(it,i);                                       % Records the Y coordinate for time of settlement for all the alive larvae
OUTDistance (i) = Rec_t_distance(it,i);                       % Records the total distance travelled for time of settlement for all the larvae
OUTgrowth (i) = Rec_growthrate(it,i);                         % Records the output of the growthrate for time of settlement for all the larvae
OUTsettled (i) = 1;                                           % Records the output whether the larvae had settled
OUTstage (i) = Rec_stage (it, i);                             % Records the stage output for time of settlement for all the larvae
elseif Rec_stage2 (it,i)== 0 && OUTgrowth(i)==0;              % Larvae that are not successful
OUTx(i) = Rec_x(it,i);                                        % Records the X coordinate for unsuccessful alive larvae
OUTy(i) = Rec_y(it,i);                                        % Records the Y coordinate for unsuccessful alive larvae
OUTiteration (i)= it;                                         % Records the iteration output for unsuccessful larvae
OUTxd(i) = Rec_x(it,i);                                       % Records the X coordinate for unsuccessful dead larvae
OUTyd(i) = Rec_y(it,i);                                       % Records the y coordinate for unsuccessful dead larvae
OUTDistance (i) = Rec_t_distance(it,i);                       % Records the total distance travelled for unsuccessful larvae
OUTgrowth (i) = Rec_growthrate(it,i);                         % Records the output growthrate for unsuccessful larvae
OUTsettled (i) = 0;                                           % Records the output for growthrate for unsuccessful larvae
OUTstage (i) = Rec_stage (it, i);                             % Records the stage for unsuccessful larvae
end                                                           % End if stage, settled and already assigned
end                                                           % End Loop for i
end                                                           % End Loop for it
Smin= 30*(OUTiteration - start);                              % Converts iteration to minutes
Shour = Smin/60;                                              % Converts minutes into hours
SDay = (Shour/24);                                            % Converts hours to days
OUTSettleday = fix(SDay);                                     % Rounds the day down i.e. 1.6 to 1 = day 1
% finds how many larvae in each patch
Celtic_Inpolgon = inpolygon(OUTx,OUTy,xpCeltic,ypCeltic);     % Output for Celtic scallop patch
C_Bay_Inpolygon = inpolygon(OUTx,OUTy,xpC_Bay,ypC_Bay);       % Output for Cardigan Bay scallop patch
Tuskar_Inpolygon = inpolygon(OUTx,OUTy,xpTuskar,ypTuskar);    % Output for Tuskar scallop patch
Trem_Bay_Inpolygon = inpolygon(OUTx,OUTy,xpTB,ypTB);          % Output for Tremadog Bay scallop patch
IOM_main_Inpolygon = inpolygon(OUTx,OUTy,xpIOM_M,ypIOM_M);    % Output for IoM main scallop patch
IOM_small_Inpolygon = inpolygon(OUTx,OUTy,xpIOM_S,ypIOM_S);   % Output for IoM top small scallop patch
N_Wales_Inpolygon = inpolygon(OUTx,OUTy,xpNW,ypNW);           % Output for North Wales scallop patch
Cornwall_Inpolygon = inpolygon(OUTx,OUTy,xpCwall,ypCwall);    % Output for Cornwall Patch scallop patch
Llyn_Inpolygon = inpolygon(OUTx,OUTy,xpLlyn,ypLlyn);          % Output for Llyn peninsula scallop patch
% output of how many larvae in each site
OUT_Celtic_Inpolgon_connect     =  sum(Celtic_Inpolgon);      % Output for Celtic scallop patch
OUT_C_Bay_Inpolygon_connect     =   sum(C_Bay_Inpolygon);     % Output for Cardigan Bay scallop patch
OUT_Tuskar_Inpolygon_connect    =   sum(Tuskar_Inpolygon);    % Output for Tuskar scallop patch
OUT_Trem_Bay_Inpolygon_connect  = sum(Trem_Bay_Inpolygon);    % Output for Tremadog Bay scallop patch
OUT_IOM_main_Inpolygon_connect  = sum(IOM_main_Inpolygon);    % Output for IoM main scallop patch
OUT_IOM_small_Inpolygon_connect = sum(IOM_small_Inpolygon);   % Output for IoM top small scallop patch
OUT_N_Wales_Inpolygon_connect   = sum(N_Wales_Inpolygon);     % Output for North Wales scallop patch
OUT_Cornwall_Inpolygon_connect  =  sum(Cornwall_Inpolygon);   % Output for Cornwall Patch scallop patch
OUT_Llyn_Inpolygon_connect      =  sum(Llyn_Inpolygon);       % Output for Llyn peninsula scallop patch

%% Movie close all
if record ==1;                                                % If to close video
close(vid);                                                   % Closes vid
end                                                           % End record
fclose('all');                                                % Closes everything

%% %%%%%%%%%%%%%%PLOT PICTURE%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clf;                                                          % Clears current figure
imageArray = uint8(2 * rand(200));                            % Basic array to set the plot up, could be made better below but works
% Make up a colormap to apply to it.
cmap(1,:) = [1 1 1];                                          % White
cmap(2,:) = [1 1 1];                                          % white
cmap(3,:) = [1 1 1];                                          % white
 imshow(imageArray, cmap);                                    % Display image with that colormap.
 set(gcf, 'Position', get(0, 'ScreenSize'));                  % Maximize figure.
clf                                                           % Clear figure so dont get overlapping figures
hold on;                                                      % Hold on allows multiple additions to same figure without over writing
whitebg('w');                                                 % White background
PTM = pcolor(long, lati,landsea');                            % Plot of the domain
set(PTM, 'EdgeColor', 'none');                                % Removes boxes
shading flat                                                  % Shading of Plot
 Celticpatch=plot(xpCeltic,ypCeltic,'k');                     % Plots Celtic scallop patch
CBaypatch=plot(xpC_Bay,ypC_Bay,'k');                          % Plots Cardigan Bay scallop patch
Tuskarpatch=plot(xpTuskar,ypTuskar,'k');                      % Plots Tuskar scallop patch
Trem_Baypatch=plot(xpTB,ypTB,'k');                            % Plots Tremadog Bay scallop patch
IOM_Mpatch=plot(xpIOM_M,ypIOM_M,'k');                         % Plots IoM main scallop patch
IOM_Spatch=plot(xpIOM_S,ypIOM_S,'k');                         % Plots IoM top small scallop patch
NWpatch=plot(xpNW,ypNW,'k');                                  % Plots North Wales scallop patch
Cornwallpatch=plot(xpCwall,ypCwall,'k');                      % Plots Cornwall Patch scallop patch
Llynpatch=plot(xpLlyn,ypLlyn,'k');                            % Plots Llyn peninsula scallop patch
plot(OUTxd, OUTyd,'xr')                                       % plots all dead larvae as red cross
plot(OUTxa, OUTya,'.k')                                       % Plots all larvae as black dots
set(gca,'color', [0 .6 0],...                                 % Changes background colour(map) to green
    'DataAspectRatio',[1 yspac/xspac 1]);                     % Plots it axis properly so not stretched        
axis([X1 X2 Y1 Y2]);                                          % Axis of full extent of the Irish sea
xlabel ('Longitude');                                         % x-axis label
ylabel('Latitude');                                           % y-axis label