% clear all;
% close all;
% clc;
%Number_of_correlated_elements=zeros(1,1545);
Max_Num=zeros(10,10);
j=1;
for i=1:10
    display(i);
    x=reshape(xi_2(i,:)',[870,870]);
    x=(((x)<93.89)&(x>0));
    x=x+0;
    x(x==0)=NaN;
    
        display(j);
        y=reshape(xi_2(j,:)',[870,870]);
        y=(((y)<93.89)&(y>0));
        y=y+0;
        y(y==0)=NaN;
        l=x==y;
        a=nnz(l);
        %Number_of_correlated_elements(i,j)=a;
        
        Max_Num(i,j)=(a/756900)*100;
        %Percentage_correlated(i,j)=((Percentage_correlated(i,j))/0.3575)*100;
    
    j=j+1;
end