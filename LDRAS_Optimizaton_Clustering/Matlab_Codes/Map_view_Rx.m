clear all;
close all;
clc;
Positions_vector_x=Positions_calc();
Positions_vector_y=Positions_vector_x;
Number_of_locations_y=numel(Positions_vector_y);
Number_of_locations_x=numel(Positions_vector_x);

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




xi_2=zeros(Number_of_locations_x,Number_of_locations_y);
a=1;
b=1;
for i=1:numel(Positions_vector_x)
    display(i);
    for j=1:numel(Positions_vector_y)
    
    
        if map(Positions_vector_x(a),Positions_vector_y(b))==0
            xi_2(i,j)=Generate_map_rx_locations(Positions_vector_x(a),Positions_vector_y(b));
            
            %xi_2=xi_2+xi_2;
        end
    b=b+1;
    end
    a=a+1;
    b=1;
end