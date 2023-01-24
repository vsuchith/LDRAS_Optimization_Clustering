function output=Path_loss_micro_build(x,y)
%Generating Map

no_of_buildings_x=10;
no_of_buildings_y=10;
width_of_buildings=140;
width_of_street=30;
res=10;
htx=1.5;
freq=800*(10^6);
MCL=0;


m=ceil((no_of_buildings_x*(width_of_buildings+width_of_street)+width_of_street)/res);
n=ceil((no_of_buildings_y*(width_of_buildings+width_of_street)+width_of_street)/res);
map=ones(m,n);
limit_map = [no_of_buildings_x*(width_of_buildings+width_of_street)+width_of_street,no_of_buildings_y*(width_of_buildings+width_of_street)+width_of_street];

%% DEFAULT PARAMETERS
hm = 1.5;
htx_eff = htx - 1.0;
hrx_eff = hm - 1.0;

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

%1) UMi: check if the transmission antenna is located in the building
if map(Location_of_transmitter_x,Location_of_transmitter_y)==0
    error('The transmission antenna must be located in the building');  
end

% Initialize output variable 
loss_outdoor = zeros(m,n,1);
loss_outdoor(:,:,1)=NaN;


%% OUTDOOR: 
% disp('Micro - Outdoor calculations');

% 
% This model ditinguishes the main street, perpedicular streets 
% and parallel streets
[j,k] = find (map == 1);

 x_bs = (Location_of_transmitter_x)*res;
 y_bs = (Location_of_transmitter_y)*res;

% Find building where the tx is located
bs_building_x_1 = (real_wall_points_x(1:2:end) <= x_bs);
bs_building_x_2 = (real_wall_points_x(2:2:end) >= x_bs);
bs_building_x = find((bs_building_x_1 & bs_building_x_2) == 1);
bs_building_y_1 = (real_wall_points_y(1:2:end) <= y_bs);
bs_building_y_2 = (real_wall_points_y(2:2:end) >= y_bs);
bs_building_y = find((bs_building_y_1 & bs_building_y_2) == 1);
 

 for l = 1:length(j)    
   
    
    x_p = (j(l))*res;
    y_p = (k(l))*res;
    
    % Find the building where the rx is located
    bs_building_x_1_p = (real_wall_points_x(1:2:end) <= x_p);
    bs_building_x_2_p = (real_wall_points_x(2:2:end) >= x_p);
    bs_building_x_p = find((bs_building_x_1_p & bs_building_x_2_p) == 1);
    bs_building_y_1_p = (real_wall_points_y(1:2:end) <= y_p);
    bs_building_y_2_p = (real_wall_points_y(2:2:end) >= y_p);
    bs_building_y_p = find((bs_building_y_1_p & bs_building_y_2_p) == 1);
 


    
    % Check if the receiver is within building
    isInbuilding = 0;
    if (~isempty(bs_building_x_p))
        if(~isempty(bs_building_y_p))
         if (bs_building_x == bs_building_x_p) && (bs_building_y == bs_building_y_p)
            isInbuilding = 1;            
         end
        end
    end
   
    
    if isInbuilding
        d_out = sqrt( (x_bs - x_p)^2 + (y_bs - y_p)^2 );
        Lout = 18.7*log10(d_out) + 46.8 + 20*log10(freq/1e9); 
        PL = Lout+20;
        loss_outdoor(j(l),k(l),1)=max(PL,MCL);
  
    end
    
       
end
% 

%% INDOOR: only for streets in LoS
%disp('Micro - Outdoor calculations');

[j,k] = find (map == 0);




for l = 1:length(j)    
   % it_current=processingState(l,length(j),it_current);    

    
    x_p = (j(l))*res;
    y_p = (k(l))*res;

    % Find the street where the rx is located
    rx_street_x_1 = (street_limits_x(1:2:end) <= x_p);
    rx_street_x_2 = (street_limits_x(2:2:end) >= x_p);
    rx_street_x = find((rx_street_x_1 & rx_street_x_2) == 1);
    rx_street_y_1 = (street_limits_y(1:2:end) <= y_p);
    rx_street_y_2 = (street_limits_y(2:2:end) >= y_p);
    rx_street_y = find((rx_street_y_1 & rx_street_y_2) == 1);
    

    if isempty(rx_street_x)
        rx_street_x=0;
    else
        rx_street_x=rx_street_x;
    end
    if isempty(rx_street_y)
        rx_street_y=0;
    else
        rx_street_y=rx_street_y;
    end
    % Check if this street is in Los of the transmitter
    isInLoS = 0;
    if (~isempty(rx_street_x))
         if (rx_street_x == bs_building_x) || (rx_street_x == bs_building_x+1)
            isInLoS = 1;            
         end
    end
    if (~isempty(rx_street_y))
         if (rx_street_y == bs_building_y) || (rx_street_y == bs_building_y+1)
            isInLoS = 1;       
         end
    end
    
    if isInLoS
        d_out = sqrt( (x_bs - x_p)^2 + (y_bs - y_p)^2 );
         Lout = 40*log10(d_out) + 7.8 - 18*log10(htx_eff) - 18*log10(hrx_eff) + 2*log10(freq/1e9); 
        PL = Lout;
        loss_outdoor(j(l),k(l),1)=max(PL,MCL);
  
    end
        
    
end
output=loss_outdoor;


            



                
















