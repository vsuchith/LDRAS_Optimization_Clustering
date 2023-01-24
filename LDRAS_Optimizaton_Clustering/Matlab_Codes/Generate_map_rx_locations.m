function loss_outdoor=Generate_map_rx_locations(x,y)
%Generating Map

no_of_buildings_x=5;
no_of_buildings_y=5;
width_of_buildings=150;
width_of_street=20;
res=1;


m=ceil((no_of_buildings_x*(width_of_buildings+width_of_street)+width_of_street)/res);
n=ceil((no_of_buildings_y*(width_of_buildings+width_of_street)+width_of_street)/res);
map=ones(m,n);
limit_map = [no_of_buildings_x*(width_of_buildings+width_of_street)+width_of_street,no_of_buildings_y*(width_of_buildings+width_of_street)+width_of_street];

%% DEFAULT PARAMETERS

%Locating Base station at the center
Location_of_transmitter_x = x;
Location_of_transmitter_y = y;

% Streets in x-direction
j = 1:m;
j = (j-1)*res;
j = mod(j,width_of_buildings+width_of_street);
k = find(j<width_of_street);
map(k,:) = 0;

% Streets in y-direction
j = 1:n;
j = (j-1)*res;
j = mod(j,width_of_buildings+width_of_street);
k = find(j<width_of_street);
map(:,k) = 0;



wp=1;
real_wall_points_x = zeros(1,no_of_buildings_x);
for r=1:no_of_buildings_x
    real_wall_points_x(wp) = (width_of_street+(r-1)*(width_of_street + width_of_buildings));
    real_wall_points_x(wp+1) = (r*(width_of_street + width_of_buildings));
    wp=wp+2;
end
wp=1;
real_wall_points_y = zeros(1,no_of_buildings_y);
for r=1:no_of_buildings_y
    real_wall_points_y(wp) = (width_of_street+(r-1)*(width_of_street + width_of_buildings));
    real_wall_points_y(wp+1) = (r*(width_of_street + width_of_buildings));
    wp=wp+2;
end

street_limits_x = [0 real_wall_points_x limit_map(1)];
street_limits_y = [0 real_wall_points_y limit_map(2)];

% 1) UMi: check if the transmission antenna is located in the street
if map(Location_of_transmitter_x,Location_of_transmitter_y)==1
    error('The transmission antenna must be located in the street');  
end

% Initialize output variable 
%loss_outdoor = zeros(m,n);
%loss_indoor = zeros(m,n);

%% OUTDOOR PS#1: 
disp('Micro - Outdoor calculations - PS#1');

[j,k] = find (map == 0);

 x_bs = (Location_of_transmitter_x)*res;
 y_bs = (Location_of_transmitter_y)*res;

for l = 1:length(j)    
   
    
    x_p = (j(l))*res;
    y_p = (k(l))*res;
    
    if (x_p==x_bs) && (y_bs==y_p)
        loss_outdoor=150;
    end
      
end

%output=reshape(loss_outdoor,[1,numel(loss_outdoor)]);
            



                
















