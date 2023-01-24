% Number of suitable clusters

%idx3 = kmeans(Percentage_correlated_d2d_tx,6,'Distance','cityblock');
idx4 = kmeans(Percentage_correlated_d2d_tx,4,'Distance','cityblock');
idx5 = kmeans(Percentage_correlated_d2d_tx,5,'Distance','cityblock');
idx6 = kmeans(Percentage_correlated_d2d_tx,6,'Distance','cityblock');
idx7 = kmeans(Percentage_correlated_d2d_tx,7,'Distance','cityblock');
idx8 = kmeans(Percentage_correlated_d2d_tx,8,'Distance','cityblock');
idx9 = kmeans(Percentage_correlated_d2d_tx,9,'Distance','cityblock');
idx10 = kmeans(Percentage_correlated_d2d_tx,10,'Distance','cityblock');
disp('Completed clustering');
[silh4,h4] = silhouette(Percentage_correlated_d2d_tx,idx4,'cityblock');
disp('Completed 1/7');
[silh5,h5] = silhouette(Percentage_correlated_d2d_tx,idx5,'cityblock');
disp('Completed 2/7');
[silh6,h6] = silhouette(Percentage_correlated_d2d_tx,idx6,'cityblock');
disp('Completed 3/7');
[silh7,h7] = silhouette(Percentage_correlated_d2d_tx,idx7,'cityblock');
disp('Completed 4/7');
[silh8,h8] = silhouette(Percentage_correlated_d2d_tx,idx8,'cityblock');
disp('Completed 5/7');
[silh9,h9] = silhouette(Percentage_correlated_d2d_tx,idx9,'cityblock');
disp('Completed 6/7');
[silh10,h10] = silhouette(Percentage_correlated_d2d_tx,idx10,'cityblock');
disp('Completed 7/7');