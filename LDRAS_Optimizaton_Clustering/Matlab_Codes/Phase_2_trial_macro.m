
clear all,
close all;
clc;

% Taken from https://www.metis2020.com/documents/simulations/ and modified
% to our requirements


% DEFAULT PARAMETERS
hm = 1.5;       % height of receiver
gammac_db=5;
R=-115;
bwidth=120;
stwidth=10;
nb_bx=7;
nb_by=7;
res=1;
hb=25;
freq=800*(10^6);
MCL=3;
% Parameters for micro-cell
no_of_buildings_x=nb_bx;
no_of_buildings_y=nb_by;
width_of_buildings=bwidth;
width_of_street=stwidth;
htx=1.5;
htx_eff = htx - 1.0;
hrx_eff = hm - 1.0;

% LAYOUT

% Size of the map: (m x n)
m = ceil((nb_bx*(bwidth+stwidth)+stwidth)/res);
n = ceil((nb_by*(bwidth+stwidth)+stwidth)/res);
map = ones(m,n);
limit_map = [nb_bx*(bwidth+stwidth)+stwidth, nb_by*(bwidth+stwidth)+stwidth];

% BS location
% if nargin~=13
%     mt = ceil(m/2);
%     nt = ceil(n/2);
% else
%     mt=round(x/res);
%     nt=round(y/res);
% end


%BS location
mt = ceil(m/2);
nt = ceil(n/2);

%Mobile tx location
Location_of_transmitter_x = 60;
Location_of_transmitter_y = 131;


% Streets in x-direction
j = 1:m;
j = (j-1)*res;
j = mod(j,bwidth+stwidth);
k = find(j<stwidth);
map(k,:) = 0;

% Streets in y-direction
j = 1:n;
j = (j-1)*res;
j = mod(j,bwidth+stwidth);
k = find(j<stwidth);
map(:,k) = 0;

% Position of building walls in x-direction and y-direction
wp=1;
real_wall_points_x = zeros(1,nb_bx);
for r=1:nb_bx
    real_wall_points_x(wp) = (stwidth+(r-1)*(stwidth + bwidth));
    real_wall_points_x(wp+1) = (r*(stwidth + bwidth));
    wp=wp+2;
end
wp=1;
real_wall_points_y = zeros(1,nb_by);
for r=1:nb_by
    real_wall_points_y(wp) = (stwidth+(r-1)*(stwidth + bwidth));
    real_wall_points_y(wp+1) = (r*(stwidth + bwidth));
    wp=wp+2;
end

street_limits_x = [0 real_wall_points_x limit_map(1)];
street_limits_y = [0 real_wall_points_y limit_map(2)];

% For Urban Macrocell -> Check if the transmission antenna is located over a building
if map(mt,nt)==0
    error('The transmission antenna must be located over a rooftop');
else
    % BS position
    x_bs = (mt-0.5)*res;
    y_bs = (nt-0.5)*res;
    
end

% 1) UMi: check if the transmission antenna is located in the street
if map(Location_of_transmitter_x,Location_of_transmitter_y)==1
    error('The transmission antenna must be located in the street');  
end

% Initialize output variable 
loss = zeros(m,n);
loss(:,:) = NaN;
xin=zeros(m,n);
xin(:,:)=NaN;


% For outdoor PS#3
disp('Macro - Outdoor calculations - PS#3');

[j,k] = find (map == 0);    % Outdoor positions
  
it_current = 0;
for l = 1:length(j)    
  
    
    % Receiver location in map [m] (center of the corresponding res x res square)
    x_p = (j(l)-0.5)*res;
    y_p = (k(l)-0.5)*res;
    
    % Distance from BS to receiver
    dist_total = sqrt((x_bs-x_p)^2+(y_bs-y_p)^2);
       
            %loss(j(l),k(l)) = 40*log10(dist_total)+7.8-18*log10(hb)-18*log10(1.5)+2*log10(freq);
               loss(j(l),k(l)) = 24.3+(35.74*log10(dist_total));
               %loss(j(l),k(l)) = 23.06+26*log(dist_total);
               xin(j(l),k(l))=loss(j(l),k(l))-gammac_db+R-14;
                
end

% For indoor PS#4
disp('Macro - Indoor calculations - PS#4');

[j,k] = find (map == 1);    % Indoor positions
it_current = 0;
for l = 1:length(j)    
    
    
    % Receiver location in map [m] (center of the corresponding res x res square)
    x_p = (j(l)-0.5)*res;
    y_p = (k(l)-0.5)*res;
    
    % Distance from BS to receiver
    dist_total = sqrt((x_bs-x_p)^2+(y_bs-y_p)^2);
       
            %loss(j(l),k(l)) = 40*log10(dist_total)+7.8-18*log10(hb)-18*log10(1.5)+2*log10(freq)+20;
            loss(j(l),k(l)) = 24.3+(35.74*log10(dist_total))+20;
            %loss(j(l),k(l)) = 23.06+26*log(dist_total);
            %xin(j(l),k(l))=loss(j(l),k(l))-gammac_db+R-14;
end

%imagesc(xin);






