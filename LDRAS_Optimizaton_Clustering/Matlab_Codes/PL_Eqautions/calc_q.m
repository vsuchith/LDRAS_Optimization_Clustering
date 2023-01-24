%function to calculate q
function q=calc_q(xin,gammad2d,efficiency)

global p_cuedb;
p_cuedb=24-30;

global p_cue;
p_cue=10^(p_cuedb/10);

global gammac_db;
gammac_db=5;

global gammac;
gammac=10^(gammac_db/10);

gammad=gammad2d;
global alp;
alp=3.5;
global R;
R=3.16e-12;
global d;
d=20;
global p_outc;
p_outc=0.05;
global p_outd;
p_outd=0.05;

r1=real((((10^(xin/10))*gammac*(1-p_outc))/(R*(exp(-gammac/R)-1+p_outc)))^(1/alp));
r2=real(((p_cue*gammad*(1-p_outd))/((10^(xin/10))*((d)^(-alp))*(exp(-gammad/((10^(xin/10))*((d)^(-alp))))-1+p_outd)))^(1/alp));
q=((800-r1)*efficiency/r2);