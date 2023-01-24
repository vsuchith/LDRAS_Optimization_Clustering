%function to calculate q
function q=calc_q_PL(xin,gammad2d,efficiency)

global p_cuedbm;
p_cuedbm=24-30;
global gammac_db;
gammac_db=5;
gammad=gammad2d;
global R;
R=-115;
global d;
d=20;
r1=(10^((gammac_db-R+xin+14-24.3)/35.74));
r2=(10^((gammad-xin+(43.7*log10(d))+p_cuedbm)/43.7));
q=((((800-r1)*efficiency)/r2));