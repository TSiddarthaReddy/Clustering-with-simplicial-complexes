function Triangle_Conductance(cluster_one,cluster_two,adj_tensor)
assoc_cluster_one_matrix = adj_tensor(cluster_one,cluster_one,cluster_one);
assoc_num_cluster_one = sum(sum(sum(assoc_cluster_one_matrix)))/6;
disp('one');disp(assoc_num_cluster_one);

assoc_cluster_two_matrix = adj_tensor(cluster_two,cluster_two,cluster_two);
assoc_num_cluster_two = sum(sum(sum(assoc_cluster_two_matrix)))/6;
disp('two');disp(assoc_num_cluster_two);

total = sum(sum(sum(adj_tensor)))/6;
disp('total');disp(total);

cut = total-assoc_num_cluster_one-assoc_num_cluster_two;
vol_1 = assoc_num_cluster_one + cut;
vol_2 = assoc_num_cluster_two + cut;

triangle_conductance = cut/min(vol_1,vol_2);
disp('TMM-conductance');disp(triangle_conductance);
end