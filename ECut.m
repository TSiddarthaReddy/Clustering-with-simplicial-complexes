function [conduc,expension,density,nassoc,ncut,minmaxcut,assoc_S,assoc_S_bar,t_conduct_2,t_nassoc_2,t_ncut_2] = ECut(G, P_no_seloop, C,C_bar,node_num,direction)

if strcmp(direction, 'directed')
    d = sum(G,2);
    % mu_v is the total number of degrees
    total_edges_degree = sum(d);
    % is C one intra-edge for one cluster?If so, temp is total number of
    % edegs within one cluster,associativity
    assoc_S = sum(sum(G(C,C)));
    % namta is vol(c_cluster),assoc+cut
    vol_S = full(sum(d(C)));
end

if strcmp(direction, 'undirected')
    tic
    d = sum(G,2);
    % mu_v is the total number of degrees
    total_edges_degree = sum(d);
    % is C one intra-edge for one cluster?If so, temp is total number of
    % edegs within one cluster,associativity
    assoc_S = sum(sum(G(C,C)));
    assoc_S_bar = sum(sum(G(C_bar,C_bar)));
    % namta is vol(c_cluster),assoc+cut
    vol_S = full(sum(d(C)));
    vol_S_bar = total_edges_degree - vol_S;
    cut = vol_S - assoc_S;
    n_c = size(C,2); 
    n_c_bar = node_num - n_c;
    public = toc;
end
tic
conduc = ECondu(cut,vol_S, vol_S_bar);
t_conduct_2 = toc;
t_conduct_2 = t_conduct_2 + public;

tic
nassoc = ENassc(assoc_S,assoc_S_bar,vol_S,vol_S_bar);
t_nassoc_2 = toc;
t_nassoc_2 = t_nassoc_2 + public;

tic
ncut = ENcut(cut,vol_S, vol_S_bar);
t_ncut_2 = toc;
t_ncut_2 = t_ncut_2 + public;

expension =  EExpension(cut,n_c,n_c_bar);
density =  EDensity(vol_S, vol_S_bar,n_c,n_c_bar);
minmaxcut = EMinmaxcut(cut,assoc_S,assoc_S_bar);
end

function conduc = ECondu(cut,vol_S, vol_S_bar)
min_vol = min(vol_S,vol_S_bar);
if min_vol == 0
    conduc = 10e6; % the graph cut returned an empty set or the whole graph.
else
    conduc = cut/min_vol;
end
end

function expension =  EExpension(cut,node_c,n_c_bar)
min_size = min(node_c,n_c_bar);
if min_size == 0 % publish is not enough
    expension = 10e6; % it may be influence relability and accuracy
else
    expension = cut/min_size;
end
end

function density =  EDensity(vol_S, vol_S_bar,node_c,n_c_bar)% TODO: the this defination of density is not heuristical
if node_c == 0 || n_c_bar == 0
    density = 0;
else
    density = max(vol_S/(node_c^2),vol_S_bar/(n_c_bar^2));
end
end

function nassoc = ENassc(assoc_S,assoc_S_bar,vol_S,vol_S_bar)
if vol_S == 0 || vol_S_bar == 0
    nassoc = 0;
else
    nassoc = assoc_S/vol_S + assoc_S_bar/vol_S_bar;
end
end

function ncut = ENcut(cut,vol_S, vol_S_bar)
if vol_S == 0 || vol_S_bar == 0
    ncut = 10e6;
else
    ncut = cut*(1/vol_S+1/vol_S_bar);
end
end
% optimal minmaxcut value is larger than 4. so it only reserves clusters
% that value is 4. That is a deviation.
function minmaxcut = EMinmaxcut(cut,assoc_S,assoc_S_bar)
if assoc_S == 0 || assoc_S_bar == 0
    minmaxcut = 10e6;
else
    minmaxcut = cut*(1/assoc_S+1/assoc_S_bar);
end
end

