% x1=-40:-10;
% r1=zeros(1,length(x1));
% for i=1:length(x1)
% %r1(i)=10^((x1(i)+6.81)/35);
% r1(i)=10^(x1(i))*log(x1(i))*log(10)+(10^(x1(i))/x1(i));
% end
% plot(r1,x1);




%f=zeros(1,length(x));

f=@(x)(((x).^(5))./((-(0.07*x)-0.04).^18));
       %((x.^(5))./(((-0.07*x)-0.04).^(18)))
q=integral(f,7,Inf);

% plot(x,f);



% x=0:0.01:7;
% f=zeros(1,length(x));
% for i=1:length(x)
%     f(i)=3.5*((x(i))^(5/18))/(-(0.07*x(i))-0.04);
% end
% plot(x,f);