function [avail_neigh, nmbr ]= avail_neigh(x,nodes)
    pos_dir1 = posdir_multi(x,length(x));
    k = ismember(pos_dir1,nodes,'rows');
    pos_dir1(k,:) = [];
    avail_neigh = pos_dir1;
    [p,l] = size(pos_dir1);
    nmbr = p;
end