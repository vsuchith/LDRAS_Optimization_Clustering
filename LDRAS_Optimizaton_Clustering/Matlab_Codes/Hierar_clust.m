function CR = Hierar_clust(S, K)

CR = [];
d = stod(S); %convert similarity matrix to distance vector
% single linkage
Z = linkage(d,'single');
CR = [CR cluster(Z,'maxclust',K)];
% complete linkage
Z = linkage(d,'complete');
CR = [CR cluster(Z,'maxclust',K)];
% average linkage
Z = linkage(d,'average');
CR = [CR cluster(Z,'maxclust',K)];