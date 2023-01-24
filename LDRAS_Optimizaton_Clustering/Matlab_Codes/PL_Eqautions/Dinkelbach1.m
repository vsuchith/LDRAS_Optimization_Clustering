clear all;
close all;
clc;
xin=-20;
amin_k=zeros(1,15);
bmin_k=zeros(1,15);
q=zeros(1,15);
fmax=zeros(1,15);
%efficiency_tmp=ones(1,15);
efficiency_tmp=[0.1523;0.2344;0.377;0.6016;0.877;1.1758;1.4766;1.9141;2.4063;2.7305;3.3223;3.9023;4.5234;5.1152;5.5547];
gammad=(10.^[-0.76182;-0.56182;-0.36222;-0.15769;0.0185;0.21667;0.40956;0.59680;0.84708;0.97837;1.24967;1.467;1.59462;1.76808;1.96808]);
for i=1:15
q(i)=calc_q(xin,gammad(i),efficiency_tmp(i));
%fmax(i)=184;
j=1;
    while j<100 
    [amin_k(i)]=vec_f_g_1(q(i),gammad(i),efficiency_tmp(i));
    fmax(i)=-bmin_k(i);
    q(i)=calc_q(amin_k(i),gammad(i),efficiency_tmp(i)); 
    j=j+1;
    end
end


    


