function[conductance] = Edge_Conductance(cluster_one,cluster_two,adj)
% only consider undirected graph(replicated edges),thus divided by two.
% calculate associativity of the first cluster
assoc_cluster_one_matrix = adj(cluster_one,cluster_one);
assoc_num_cluster_one = sum(sum(assoc_cluster_one_matrix))/2;

% calculate associativity of the second cluster
assoc_cluster_two_matrix = (adj(cluster_two,cluster_two));
assoc_num_cluster_two = sum(sum(assoc_cluster_two_matrix))/2;

% no replicated for cut,thus no divide by 2;
cut_cluster_matrix = adj(cluster_one,cluster_two);
cut_num = sum(sum(cut_cluster_matrix));

% vol
vol_1 = assoc_num_cluster_one + cut_num;
vol_2 = assoc_num_cluster_two + cut_num;

total = sum(sum(adj));

conductance = cut_num/min(vol_1,vol_2);
disp('conductance');disp(conductance)
end

