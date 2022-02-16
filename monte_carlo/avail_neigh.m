function [avail_neigh, nmbr ]= avail_neigh(x,nodes)
pos_dir1 = posdir_2d(x,1);
k = ismember(pos_dir1,nodes,'rows');
pos_dir1(k,:) = [];
avail_neigh = pos_dir1;
[p,l] = size(pos_dir1);
nmbr = p;
end