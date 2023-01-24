function f_c = new_optobj(xin,gammac_db,R,p_cuedb,gammad,d)

r1=(10^((gammac_db-R+xin+14-24.3)/35.74));
r2=(10^((gammad-xin+(43.7*log10(d))+p_cuedb)/43.7));
% r1=real((((10^(xin(1)/10))*gammac*(1-xin(2)))/(R*(exp(-gammac/R)-1+xin(2))))^(1/alp));
% r2=real(((p_cue*gammad*(1-xin(3)))/((10^(xin(1)/10))*((d)^(-alp))*(exp(-gammad/((10^(xin(1)/10))*((d)^(-alp))))-1+xin(3))))^(1/alp));
f_c=(r1+r2);
end