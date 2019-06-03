function [conf_mat, acc, bel_new] = accuracy(label, bel, m)

% This functions tries out all the possibles permutations of the labels
% 'bel' in order to much the true labels 'label' in the confusion matrix.
% label: Nx1 matrix, true labels
% bel:   Nx1 matrix, clustering
% m: number of clusters

P = perms(1:m);
[num_perm, ~] = size(P);
conf_mat = confusionmat(label, bel);

p_idx = 0; % index of the correct permutation
S = 0;
for j=1:num_perm
   D = sum(diag(conf_mat(:,P(j,:))));
   
   if (D>S)
       S=D;
       p_idx=j;
   end

end

bel_new = changem(bel, P(p_idx, :)', (1:m));
conf_mat = conf_mat(:,P(p_idx,:));
acc = S / sum(sum(conf_mat));
end