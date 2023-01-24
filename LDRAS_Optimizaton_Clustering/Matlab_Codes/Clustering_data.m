%to mask

Loss_Macro=Phase_2_trial_macro();
Loss_Macro=Loss_Macro<103;
indices=find(Loss_Macro);


clustering_d1=zeros(350,29929);
clustering_c1=zeros(350,29929);
Path_loss_2d=reshape(permute(Path_loss,[3,2,1]),[],29929);


Path_loss_pa_vue=Path_loss_2d<93.89;
Path_loss_pa_vue=double(Path_loss_pa_vue);
for j=1:length(indices)
    Path_loss_pa_vue(indices(j),:)=NaN;
end


Path_loss_pa_cue=Path_loss_2d>104.58;
nan_log_mat=isnan(Path_loss_2d);
nan_double_mat=double(nan_log_mat);

Path_loss_pa_cue=double(Path_loss_pa_cue);
Path_loss_pa_cue=Path_loss_pa_cue+nan_double_mat;
for j=1:length(indices)
    Path_loss_pa_cue(indices(j),:)=NaN;
end


rng(1);

 parfor i=1:350
 display(i);
 num_of_clusters=randi([36,40],1);
 clustering_c1(i,:)=kmeans(Path_loss_pa_cue,num_of_clusters,'Distance','correlation','EmptyAction','singleton','Start','sample');
 clustering_d1(i,:)=kmeans(Path_loss_pa_vue,num_of_clusters,'Distance','correlation','Start','sample','EmptyAction','singleton');
end
save('clustering_c1_old.mat','clustering_c1','-v7.3');
save('clustering_d1_old.mat','clustering_d1','-v7.3');
