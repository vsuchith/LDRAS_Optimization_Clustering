%Generating Map
clear all;
close all;
clc;
no_of_buildings_x=7;
no_of_buildings_y=7;
width_of_buildings=120;
width_of_street=20;
res=1;
htx=1.5;
freq=800*(10^6);
MCL=3;
hs=5;
m=ceil((no_of_buildings_x*(width_of_buildings+width_of_street)+width_of_street)/res);
n=ceil((no_of_buildings_y*(width_of_buildings+width_of_street)+width_of_street)/res);
map=ones(m,n);
limit_map = [no_of_buildings_x*(width_of_buildings+width_of_street)+width_of_street,no_of_buildings_y*(width_of_buildings+width_of_street)+width_of_street];

%% DEFAULT PARAMETERS
hm = 1.5;
htx_eff = htx - 1.0;
hrx_eff = hm - 1.0;
gammad_db=7;
p_cue=5;
%Locating Base station at the center
Location_of_transmitter_x = 2;
Location_of_transmitter_y = 1;

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
loss = zeros(m,n);
loss(:,:) = NaN;

%% OUTDOOR PS#1: 
disp('Micro - Outdoor calculations - PS#1');

% The index of each vector are ordered in the following directions:
% 1 = x-positive (East)
% 2 = y-positive (North)
% 3 = x-negative (West)
% 4 = y-negative (South)

% This model ditinguishes the main street, perpedicular streets 
% and parallel streets
[j,k] = find (map == 0);

 x_bs = (Location_of_transmitter_x - 0.5)*res;
 y_bs = (Location_of_transmitter_y - 0.5)*res;

% Find direction and center of the street where the BS is located
bs_street_x_1 = (street_limits_x(1:2:end) <= x_bs);
bs_street_x_2 = (street_limits_x(2:2:end) >= x_bs);
bs_street_x = find((bs_street_x_1 & bs_street_x_2) == 1);
bs_street_y_1 = (street_limits_y(1:2:end) <= y_bs);
bs_street_y_2 = (street_limits_y(2:2:end) >= y_bs);
bs_street_y = find((bs_street_y_1 & bs_street_y_2) == 1);

if ~isempty(bs_street_x)
    if ~isempty(bs_street_y)
        % The transmitter is located in a intersection of perpendicular
        % streets
        bs_str = 2;
        bs_str_centre(1) = (street_limits_x(2*bs_street_x) + street_limits_x(2*bs_street_x-1))/2;
        bs_str_centre(2) = (street_limits_y(2*bs_street_y) + street_limits_y(2*bs_street_y-1))/2;        
    else
        % The transmitter is located in a street in y-direction
        bs_str = 1;
        bs_str_centre(1) = (street_limits_x(2*bs_street_x) + street_limits_x(2*bs_street_x-1))/2;
        bs_str_centre(2) = -1;        
    end
else
    % The transmitter is located in a street in x-direction
    bs_str = 0;
    bs_str_centre(1) = -1;
    bs_str_centre(2) = (street_limits_y(2*bs_street_y) + street_limits_y(2*bs_street_y-1))/2;        
end
d1=zeros(m,n);
d1_1=zeros(m,n);
d1_2=zeros(m,n);
d2=zeros(m,n);
d2_1=zeros(m,n);
d2_2=zeros(m,n);

it_current = 0;
for l = 1:length(j)    
    it_current=processingState(l,length(j),it_current);
    
    x_p = (j(l)-0.5)*res;
    y_p = (k(l)-0.5)*res;
    
    % Find the street where the rx is located
    rx_street_x_1 = (street_limits_x(1:2:end) <= x_p);
    rx_street_x_2 = (street_limits_x(2:2:end) >= x_p);
    rx_street_x = find((rx_street_x_1 & rx_street_x_2) == 1);
    rx_street_y_1 = (street_limits_y(1:2:end) <= y_p);
    rx_street_y_2 = (street_limits_y(2:2:end) >= y_p);
    rx_street_y = find((rx_street_y_1 & rx_street_y_2) == 1);
   % rx_str_centre=zeros(1,2);
    if ~isempty(rx_street_x)
        if ~isempty(rx_street_y)
            % The receiver is located in perpendicular street            
            rx_str = 2;  
            rx_str_centre(1) = (street_limits_x(2*rx_street_x) + street_limits_x(2*rx_street_x-1))/2;
            rx_str_centre(2) = (street_limits_y(2*rx_street_y) + street_limits_y(2*rx_street_y-1))/2;        
        else
            % The receiver is located in a street in y-direction
            rx_str = 1;
            rx_str_centre(1) = (street_limits_x(2*rx_street_x) + street_limits_x(2*rx_street_x-1))/2;
            rx_str_centre(2) = -1;        
        end
    else
        % The receiver is located in a street in y-direction
        rx_str = 0;
        rx_str_centre(1) = -1;
        if isempty(rx_street_y)
            rx_str_centre(2)=0;
        else
        rx_str_centre(2) =(street_limits_y(2*rx_street_y) + street_limits_y(2*rx_street_y-1))/2;        
        end
    end


% Rect from BS to point
    if (x_p - x_bs) ~= 0
        slope = (y_p - y_bs)/(x_p - x_bs);
        ordi = (x_p*y_bs - x_bs*y_p)/(x_p-x_bs);
        normal_to_x = 0;
    else
        % Rect perpendicular to x-axis
        normal_to_x = 1;        
    end
    
    % Check if the receiver is in LoS
    isInLoS = isLoS(j(l),k(l),Location_of_transmitter_x ,Location_of_transmitter_y,x_p,y_p,x_bs,y_bs,real_wall_points_x,real_wall_points_y,slope,ordi,normal_to_x);
    %isInLoS = isLoS(j(l),k(l),mt,nt,x_p,y_p,x_bs,y_bs,real_wall_points_x,real_wall_points_y,slope,ordi,normal_to_x);
    
    if isInLoS
        % Receiver in the main street (LoS)
        d1(j(l),k(l)) = sqrt((x_p-x_bs)^2 + (y_p-y_bs)^2);
        Llos = 40*log10(d1(j(l),k(l))) + 7.8 - 18*log10(htx_eff) - 18*log10(hrx_eff) + 2*log10(freq/1e9);
        loss(j(l),k(l)) = max(Llos,MCL);  
        
    else
        % Check if the receiver is either in a perpendicular or a paralel
        % street
        parallel = false;
        if ((rx_str == bs_str) && (rx_str < 2))
            parallel = true;      
        end
        
        if parallel
            loss(j(l),k(l)) = NaN;
        else            
            if rx_str * bs_str < 4
                switch bs_str
                    case 0
                        d1(j(l),k(l)) = abs(x_bs-rx_str_centre(1));
                        d2(j(l),k(l)) = abs(y_p-bs_str_centre(2));
                    case 1
                        d1(j(l),k(l)) = abs(y_bs-rx_str_centre(2));
                        d2(j(l),k(l)) = abs(x_p-bs_str_centre(1));
                    case 2
                        if rx_str == 0
                            d1(j(l),k(l)) = abs(y_bs-rx_str_centre(2));
                            d2(j(l),k(l)) = abs(x_p-bs_str_centre(1));
                        else
                            d1(j(l),k(l)) = abs(x_bs-rx_str_centre(1));
                            d2(j(l),k(l)) = abs(y_p-bs_str_centre(2));
                        end
                end
                Llos1 = 40*log10(d1(j(l),k(l))) + 7.8 - 18*log10(htx_eff) - 18*log10(hrx_eff) + 2*log10(freq/1e9);
                Llos2 = 40*log10(d2(j(l),k(l))) + 7.8 - 18*log10(hrx_eff) - 18*log10(hrx_eff) + 2*log10(freq/1e9);
                nj1 = max(2.8-0.0024*d1(j(l),k(l)),1.84);
                nj2 = max(2.8-0.0024*d2(j(l),k(l)),1.84);
                PL1 = Llos1 + 17.9 - 12.5*nj1 + 10*nj1*log10(d2(j(l),k(l))) + 3*log10(freq/1e9);
                PL2 = Llos2 + 17.9 - 12.5*nj2 + 10*nj2*log10(d1(j(l),k(l))) + 3*log10(freq/1e9);
                PL = min(PL1,PL2);                  
            else
                d1_1(j(l),k(l)) = abs(x_bs-rx_str_centre(1));
                d2_1(j(l),k(l)) = abs(y_p-bs_str_centre(2));
                Llos1_1 = 40*log10(d1_1(j(l),k(l))) + 7.8 - 18*log10(htx_eff) - 18*log10(hrx_eff) + 2*log10(freq/1e9);
                Llos2_1 = 40*log10(d2_1(j(l),k(l))) + 7.8 - 18*log10(hrx_eff) - 18*log10(hrx_eff) + 2*log10(freq/1e9);
                nj1_1 = max(2.8-0.0024*d1_1(j(l),k(l)),1.84);
                nj2_1 = max(2.8-0.0024*d2_1(j(l),k(l)),1.84);
                PL1_1 = Llos1_1 + 17.9 - 12.5*nj1_1 + 10*nj1_1*log10(d2_1(j(l),k(l))) + 3*log10(freq/1e9);
                PL2_1 = Llos2_1 + 17.9 - 12.5*nj2_1 + 10*nj2_1*log10(d1_1(j(l),k(l))) + 3*log10(freq/1e9);
                PL_1 = min(PL1_1,PL2_1);

                d1_2(j(l),k(l)) = abs(y_bs-rx_str_centre(2));
                d2_2(j(l),k(l)) = abs(x_p-bs_str_centre(1));
                Llos1_2 = 40*log10(d1_2(j(l),k(l))) + 7.8 - 18*log10(htx_eff) - 18*log10(hrx_eff) + 2*log10(freq/1e9);
                Llos2_2 = 40*log10(d2_2(j(l),k(l))) + 7.8 - 18*log10(hrx_eff) - 18*log10(hrx_eff) + 2*log10(freq/1e9);
                nj1_2 = max(2.8-0.0024*d1_2(j(l),k(l)),1.84);
                nj2_2 = max(2.8-0.0024*d2_2(j(l),k(l)),1.84);
                PL1_2 = Llos1_2 + 17.9 - 12.5*nj1_2 + 10*nj1_2*log10(d1_2(j(l),k(l))) + 3*log10(freq/1e9);
                PL2_2 = Llos2_2 + 17.9 - 12.5*nj2_2 + 10*nj2_2*log10(d2_2(j(l),k(l))) + 3*log10(freq/1e9);
                PL_2 = min(PL1_2,PL2_2);

                PL = min(PL_1,PL_2);
            end
            loss(j(l),k(l)) = max(PL,MCL);  
           
        end
    end    
end


%% INDOOR: only for buiding in LoS
disp('Micro - Indoor calculations - PS#2');

% The index of each vector are ordered in the following directions:
% 1 = x-positive (East)
% 2 = y-positive (North)
% 3 = x-negative (West)
% 4 = y-negative (South)
loss_indoor= zeros(m,n);
%loss_indoor(:,:) = NaN;
[j,k] = find (map == 1);

normal_vect = [1 0; 0 1; -1 0; 0 -1];

it_current = 0;
for l = 1:length(j)    
    it_current=processingState(l,length(j),it_current);    

    
    x_p = (j(l) - 0.5)*res;
    y_p = (k(l) - 0.5)*res;

    % Find the build where the rx is located
    rx_build_x_1 = (real_wall_points_x(1:2:end) <= x_p);
    rx_build_x_2 = (real_wall_points_x(2:2:end) >= x_p);
    rx_build_x = find((rx_build_x_1 & rx_build_x_2) == 1);
    rx_build_y_1 = (real_wall_points_y(1:2:end) <= y_p);
    rx_build_y_2 = (real_wall_points_y(2:2:end) >= y_p);
    rx_build_y = find((rx_build_y_1 & rx_build_y_2) == 1);
    

    if isempty(rx_build_x)
        rx_build_x=0;
    else
        rx_build_x=rx_build_x;
    end
    if isempty(rx_build_y)
        rx_build_y=0;
    else
        rx_build_y=rx_build_y;
    end
    % Check if this building is in the street of the transmitter
    isInLoS = 0;
    wall_sighted = zeros(4,1);
    wall_points = -1*ones(4,2);
    d_to_wall = -1*ones(4,1);
    if (~isempty(bs_street_x))
        if (rx_build_x == bs_street_x) || (rx_build_x == bs_street_x - 1)|| (rx_build_x == bs_street_x - 2)|| (rx_build_x == bs_street_x - 3)
            isInLoS = 1;            
            if (rx_build_x + 1 == bs_street_x)
                wall_sighted(1) = 1;
                wall_points(1,:) = [real_wall_points_x(2*rx_build_x), y_p];
                d_to_wall(1) = abs(x_p - real_wall_points_x(2*rx_build_x));
            else
                wall_sighted(3) = 1;
                wall_points(3,:) = [real_wall_points_x(2*rx_build_x - 1), y_p];
                d_to_wall(3) = abs(x_p - real_wall_points_x(2*rx_build_x - 1));
            end
        end        
    end
    if (~isempty(bs_street_y))
        if (rx_build_y == bs_street_y) || (rx_build_y == bs_street_y - 1)|| (rx_build_y == bs_street_y - 2)|| (rx_build_y == bs_street_y - 3)
            isInLoS = 1;            
            if (rx_build_y + 1 == bs_street_y)
                wall_sighted(2) = 1;
                wall_points(2,:) = [x_p real_wall_points_y(2*rx_build_y)];
                d_to_wall(2) = abs(y_p - real_wall_points_y(2*rx_build_y));
            else
                wall_sighted(4) = 1;
                wall_points(4,:) = [x_p real_wall_points_y(2*rx_build_y - 1)];
                d_to_wall(4) = abs(y_p - real_wall_points_y(2*rx_build_y - 1));
            end
        end       
    end
    
    if isInLoS
        if sum(wall_sighted) < 2
            num_opt = 1;
        else
            % Two options because the transmitter is in a intersection
            num_opt = 2;
        end
        
        if num_opt == 1
            idx = find(wall_sighted);
            d_in = d_to_wall(idx);           
            d_out = sqrt( (x_bs - wall_points(idx,1))^2 + (y_bs - wall_points(idx,2))^2 );
            
            
            
            Lout = 40*log10(d_out + d_in) + 7.8 - 18*log10(htx_eff) - 18*log10(hrx_eff) + 2*log10(freq/1e9); 
            Lout = max(Lout,MCL);      
            PL = Lout;
            loss(j(l),k(l))=max(PL,MCL)+20;
        else
            
            % First option
            idx = find(wall_sighted,1,'first');
            
            d_in1 = d_to_wall(idx);           
            d_out1 = sqrt( (x_bs - wall_points(idx,1))^2 + (y_bs - wall_points(idx,2))^2 );
            
            % Second option
            idx = find(wall_sighted,1,'last');
            
            d_in2 = d_to_wall(idx);           
            d_out2 = sqrt( (x_bs - wall_points(idx,1))^2 + (y_bs - wall_points(idx,2))^2 );
            
            Lout1 = 40*log10(d_out1 + d_in1) + 7.8 - 18*log10(htx_eff) - 18*log10(hrx_eff) + 2*log10(freq/1e9); 
            Lout1 = max(Lout1,MCL);               
            PL1 = Lout1;
                    
            Lout2 = 40*log10(d_out2 + d_in2) + 7.8 - 18*log10(htx_eff) - 18*log10(hrx_eff) + 2*log10(freq/1e9); 
            Lout2 = max(Lout2,MCL);
            PL2 = Lout2;                    
                    
            PL = min(PL1,PL2);
                    
            loss(j(l),k(l))=max(PL,MCL)+20;
        end
   
    end
    
end

% imagesc(loss);
% C=colormap;
% L=size(C,1);
% loss_scale=round(interp1(linspace(min(loss(:)),max(loss(:)),L),1:L,loss));
% H = reshape(C(loss_scale,:),[size(loss_scale) 3]); 
% image(H);

%hist(loss);
% loss_1=rgb2gray(loss);
% hist(loss_1);
                
















