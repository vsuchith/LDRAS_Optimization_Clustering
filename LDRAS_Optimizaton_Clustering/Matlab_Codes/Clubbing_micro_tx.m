clear all;
close all;
clc;


%Generating Map

no_of_building s_x=5;
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

%% to calculate the percentage correlated

% Number_of_locations_y=2;
% Number_of_locations_x=2;
% Number_of_correlated_elements=zeros(Number_of_locations_x,Number_of_locations_y);
% Percentage_correlated=zeros(Number_of_locations_x,Number_of_locations_y);

%xi=zeros(Number_of_locations,739600);
% d_x=2;
% d_y=1;
% for i=1:Number_of_locations_x
%     display(i);
%     for j=1:Number_of_locations_y
%     
%     
%     if map(d_x,d_y)==0
%     xi_2=Los_indoor(d_x,d_y);
%     coco=reshape(xi_2,[870,870]);
%     l=co==coco;
%     a=nnz(l);
%     Number_of_correlated_elements(i,j)=a;
%     Percentage_correlated(i,j)=(Number_of_correlated_elements(i,j)/756900)*100;
%     
%     end
%     d_y=d_y+10;
%     end
%     d_x=d_x+10;
%     d_y=1;
% end

%% To calculate the path loss for each rx location


Positions_vector_x=Positions_calc();
Positions_vector_y=Positions_vector_x;
Number_of_locations_y=numel(Positions_vector_y);
Number_of_locations_x=numel(Positions_vector_x);

%xi_2=zeros(Number_of_locations_x,Number_of_locations_y);
xi_2=zeros(618,756900);
a=4;
b=1;
c=1;
for i=1:6
    display(i);
    for j=1:numel(Positions_vector_y)
        display(j);
    
        if map(Positions_vector_x(a),Positions_vector_y(b))==0
            xi_2(c,:)=Los_indoor(Positions_vector_x(a),Positions_vector_y(b));
            c=c+1;
            
        end
    
    b=b+1;
    end
    a=a+1;
    b=1;
end
%save('first_1-4','xi_2');

%% miscellaneous
% for i=2:numel(Number_of_correlated_elements)
%     co=reshape(xi(1,:),[860,860]);
%     coco=reshape(xi(i,:),[860,860]);
%     l=co==coco;
%     a=nnz(l);
%     Number_of_correlated_elements(1,i)=a;
%     Percentage_correlated(1,i)=(Number_of_correlated_elements(1,i)/739600)*100;
% end

%xin_2=Los_indoor(20,1);
%xin_2=Los_indoor(739,1);
% xin_3=Los_indoor(490,350);
% xin_4=Los_indoor(490,1);
% xin_5=Los_indoor(650,60);
%xi=xin_2+xin_1+xin_3+xin_4+xin_5;


% idx=kmeans(xi,20);
% imagesc(xi);
% colormap('gray')
% C=colormap;
% L=size(C,1);
% loss_scale=round(interp1(linspace(min(xi(:)),max(xi(:)),L),1:L,xi));
% H = reshape(C(loss_scale,:),[size(loss_scale) 3]); 
% z=image(H);
% [lev]=MET(z);

% 
% if(xin_1==xin_2)
%     v=0;
% else
%     v=1;
% end