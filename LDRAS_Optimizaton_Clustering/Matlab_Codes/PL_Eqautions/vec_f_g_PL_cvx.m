function [y,r1,r2] = vec_f_g_PL_cvx(q1,gammad2d,efficiency)

% r1_vec=;
% r2_vec=zeros(1,15);
q_temp=q1;
gammad=gammad2d;
gammac_db=5;
d=20;
R=-115;
p_cuedb=24-30;
% f_c=@(xin)optobj_PL(xin,q_temp,gammac_db,R,p_cuedb,gammad,d,efficiency_temp);
% options=optimoptions('fmincon','Hessian','bfgs');
% [a,b]=fmincon(f_c,-29,[],[],[],[],-40,-26,[],options);
cvx_begin
    variable xin
    r1=(10^((gammac_db-R+xin+14-24.3)/35.74));
    r2=(10^((gammad-xin+(43.7*log10(d))+p_cuedb)/43.7));
    f_c=-(((800-r1)*efficiency-(q_temp*r2)));
    minimize(f_c)
    subject to
        -40 <= xin <= -10
        (10^((gammac_db-R+xin+14-24.3)/35.74))<=200
        (10^((gammad-xin+(43.7*log10(d))+p_cuedb)/43.7))<=216
        
cvx_end
y=xin;
r1=(10^((gammac_db-R+y+14-24.3)/35.74));
r2=(10^((gammad-y+(43.7*log10(d))+p_cuedb)/43.7));

end


