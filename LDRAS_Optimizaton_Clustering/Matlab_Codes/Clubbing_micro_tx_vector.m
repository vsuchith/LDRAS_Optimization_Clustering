clear all;
close all;
clc;
Number_of_locations=30;
Number_of_correlated_elements=zeros(1,Number_of_locations);
Percentage_correlated=zeros(1,Number_of_locations);
xi=zeros(Number_of_locations,739600);
d=2;
for i=1:Number_of_locations
    display(i);
xi(i,:)=Los_indoor(2,d);
d=d+15;
end
for i=2:numel(Number_of_correlated_elements)
    co=reshape(xi(1,:),[860,860]);
    coco=reshape(xi(i,:),[860,860]);
    l=co==coco;
    a=nnz(l);
    Number_of_correlated_elements(1,i)=a;
    Percentage_correlated(1,i)=(Number_of_correlated_elements(1,i)/739600)*100;
end

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