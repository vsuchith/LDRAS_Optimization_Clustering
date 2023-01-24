
%to replace NaNs in final clustered information
clust=[clustering_c1;clustering_d1];
b=max(clust,[],2);

for i=1:size(clust,1)

[~,y1]=find(isnan(clust(i,:)));

clust(i,y1(:))=b(i)+1;

end
