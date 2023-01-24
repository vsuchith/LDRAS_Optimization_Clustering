function [a,b] = vec_f_g_PL(q1,gammad2d,efficiency)

% r1_vec=zeros(1,15);
% r2_vec=zeros(1,15);
q_temp=q1;
efficiency_temp=efficiency;
gammad=gammad2d;
gammac_db=5;
d=20;
R=-115;
p_cuedb=24-30;
f_c=@(xin)optobj_PL(xin,q_temp,gammac_db,R,p_cuedb,gammad,d,efficiency_temp);
options=optimoptions('fmincon','Hessian','bfgs');
[a,b]=fmincon(f_c,-29,[],[],[],[],-40,-26,[],options);


end


