b=max(clustering_c1,[],2);
l=find(b==42);
x=find(b==43);
y=find(b==44);
z=find(b==45);
all_ind=[l;x;y;z];
for i=1:350
clustering_c1(all_ind,:)=0;
end
clustering_c1(~any( clustering_c1, 2 ), : ) = [];
%append