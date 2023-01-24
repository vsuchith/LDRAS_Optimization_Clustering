
function output=Path_loss_macro()


% DEFAULT PARAMETERS
hm = 1.5;       % height of receiver

bwidth=140;
stwidth=30;
nb_bx=10;
nb_by=10;
res=10;
freq=800*(10^6);
MCL=0;
% Parameters for micro-cell
no_of_buildings_x=nb_bx;
no_of_buildings_y=nb_by;
width_of_buildings=bwidth;
width_of_street=stwidth;
htx=1.5;
htx_eff = htx - 1.0;
hrx_eff = hm - 1.0;



% Size of the map: (m x n)
m = ceil((nb_bx*(bwidth+stwidth)+stwidth)/res);
n = ceil((nb_by*(bwidth+stwidth)+stwidth)/res);
map = ones(m,n);
limit_map = [nb_bx*(bwidth+stwidth)+stwidth, nb_by*(bwidth+stwidth)+stwidth];


%BS location
mt = 80;
nt = 92;


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


if map(mt,nt)==0
    error('The transmission antenna must be located over a rooftop');
else
    % BS position
    x_bs = (mt)*res;
    y_bs = (nt)*res;
    
end


% Initialize output variable 
loss = zeros(m,n);
loss(:,:) = NaN;



% For outdoor PS#3
disp('Macro - Outdoor calculations');

[j,k] = find (map == 0);    % Outdoor positions

for l = 1:length(j)    
  
    
    % Receiver location in map [m] (center of the corresponding res x res square)
    x_p = (j(l))*res;
    y_p = (k(l))*res;
    
    % Distance from BS to receiver
    dist_total = sqrt((x_bs-x_p)^2+(y_bs-y_p)^2);
    loss(j(l),k(l)) = 24.3+(35.74*log10(dist_total));
              
                
end

% For indoor PS#4
disp('Macro - Indoor calculations ');

[j,k] = find (map == 1);    % Indoor positions

for l = 1:length(j)    
    
    
    % Receiver location in map [m] (center of the corresponding res x res square)
    x_p = (j(l))*res;
    y_p = (k(l))*res;
    
    % Distance from BS to receiver
    dist_total = sqrt((x_bs-x_p)^2+(y_bs-y_p)^2);
    loss(j(l),k(l)) = 24.3+(35.74*log10(dist_total));
           
end
output=loss;
%imagesc(xin);






