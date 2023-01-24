clear all;
close all;
clc;
xin=-20;
q=zeros(1,100);
%calculate q with present xin
q(1)=calc_q(xin,2.5,1.4);
%call function to solve maximisation problem
amin_k=zeros(1,100);
bmin_k=zeros(1,100);

for i=1:100
[amin_k(i),bmin_k(i)]=vec_f_g_1(q(i),2.5,1.4);
qmax=-bmin_k;
if i==1
q(i)=calc_q(amin_k(i),2.5,1.4);
else
q(i+1)=calc_q(amin_k(i),2.5,1.4);   
end
end

    


