clear all;
close all;
clc;

%Generating Map

no_of_buildings_x=10;
no_of_buildings_y=10;
width_of_buildings=140;
width_of_street=30;
res=10;
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

%% To calculate the path loss for each rx location

Size_of_Map=size(map,2);
Positions_vector_x=(1:Size_of_Map);
Positions_vector_y=(1:Size_of_Map);
Path_loss=zeros(Size_of_Map,Size_of_Map,Size_of_Map*(Size_of_Map));
a=1;
b=1;
c=1;
for i=1:numel(Positions_vector_x)
    display(i);
    for j=1:numel(Positions_vector_y)
        
    
            if map(Positions_vector_x(a),Positions_vector_y(b))==0
            Path_loss(:,:,c)=Los_indoor(Positions_vector_x(a),Positions_vector_y(b));
            
            
            elseif map(Positions_vector_x(a),Positions_vector_y(b))==1
              Path_loss(:,:,c)=Trial_out_ind(Positions_vector_x(a),Positions_vector_y(b));
            end
        
    c=c+1;
    b=b+1;
    end
    a=a+1;
    b=1;
end



