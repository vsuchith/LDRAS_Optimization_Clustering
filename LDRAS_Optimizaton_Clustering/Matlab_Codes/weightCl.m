function wcl = weightCl(E) 


N = size(E,1); %no of data points
no_allcl = max(max(E));
pc = zeros(N,no_allcl); % matrix indicates if data point belongs to the cluster (1=y, 0=n), row=data, col = cluster
for i=1:N
    pc(i,E(i,:))=1; % pc(i,j) = 1 if data i belongs to cluster j
end
display('To have an idea of the time it takes to simulate')
%find number of shared data points for each pair of clusters ==> intersect/union
wcl = zeros(no_allcl,no_allcl);
for i=1:no_allcl-1
    
    todisp=[num2str(i),'outof',num2str(no_allcl-1)];
    display(todisp);
    for ii=i+1:no_allcl
        tmp = sum((pc(:,i)+pc(:,ii))>0);
        if tmp > 0
            wcl(i,ii) = sum((pc(:,i)+pc(:,ii))==2) / tmp; %intersection/union
        end
    end
end

display('Time taking part is over; near to the result')
wcl = wcl + wcl';