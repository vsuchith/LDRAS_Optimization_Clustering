% clear all;
% close all;
% clc;
%Number_of_correlated_elements=zeros(10609,10609);
Percentage_correlated_d2d_tx=zeros(1,1000);
Percentage_correlated_cue=zeros(1,1000);
for i=1:1
    display(i);
    x=reshape(xi_2(i,:)',[870,870]);
    x=(((x)<93.89)&(x>0));
    x_cue=((x>104.58));
    %x=x+0;
    %x(x==0)=NaN;
    for j=1:1000
        display(j);
        y=reshape(xi_2(j,:)',[870,870]);
        y=(((y)<93.89)&(y>0));
        y_cue=((y>104.58));
        %y=y+0;
        %y(y==0)=NaN;
        l=x==y;
        m=x_cue==y_cue;
        a=nnz(l);
        b=nnz(m);
        %Number_of_correlated_elements(i,j)=a;
        
        Percentage_correlated_d2d_tx(i,j)=(a/756900)*100;
        Percentage_correlated_cue(i,j)=(b/756900)*100;
        %Percentage_correlated(i,j)=((Percentage_correlated(i,j))/0.3575)*100;
    end
    j=1;
end