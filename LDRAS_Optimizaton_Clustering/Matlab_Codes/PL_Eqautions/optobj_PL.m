function f_c = optobj_PL(xin,q_temp,gammac,R,p_cue,gammad,d,efficiency)

r1=(10^((gammac-R+xin+14-24.3)/35.74));
r2=(10^((gammad-xin+(43.7*log10(d))+p_cue)/43.7));
f_c=-(((800-r1)*efficiency-(q_temp*r2)));

end