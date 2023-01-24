clear all;
clc;
handle=@vec_f_g;
%x_o=-6;
%a=fmincon(handle,x_o,[],[],[],[],-10,20);
options=optimoptions('fmincon','Hessian','lbfgs');
%a=minConf_TMP(handle,[-6 25],[-10 20],[20 30]);
[a,b]=fmincon(handle,-6,[],[],[],[],-10,20,[],options);


%[a,b,c]=projbfgs(x_o,handle,UB,LB);