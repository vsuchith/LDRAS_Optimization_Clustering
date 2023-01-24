function final_vector=Positions_calc()
positions_between=ones(5,19);
vector=[2 10 18 172 180 188 342 350 358 512 520 528 682 690 698 852 860 868];
a=1;
j=1;
for i=1:numel(vector)
    if  (mod(i,3)==0) && (a<=5)
        positions_between(j,:)=linspace(vector(i),vector(i+1),19);
        a=a+1;
        j=j+1;
    end
end
final_vector=[2, 10, 18, positions_between(1,(2:end-1)),172, 180, 188,positions_between(2,(2:end-1)), 342, 350, 358,positions_between(3,(2:end-1)), 512, 520, 528,positions_between(4,(2:end-1)), 682, 690, 698,positions_between(5,(2:end-1)), 852, 860, 868];
final_vector=ceil(final_vector);