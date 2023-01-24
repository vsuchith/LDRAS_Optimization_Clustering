function f_c = optobj(xin,q_temp,gammac,p_outc,R,alp,p_cue,gammad,p_outd,d,efficiency)


r1=real((((10^(xin/10))*gammac*(1-p_outc))/(R*(exp(-gammac/R)-1+p_outc)))^(1/alp));
r2=real(((p_cue*gammad*(1-p_outd))/((10^(xin/10))*((d)^(-alp))*(exp(-gammad/((10^(xin/10))*((d)^(-alp))))-1+p_outd)))^(1/alp));
%f_c=-((800-r1).*(1/r2))*efficiency;
f_c=-(((800-r1)*efficiency-(q_temp*r2)));
%f_c=(r1+5*r2)/efficiency;
end