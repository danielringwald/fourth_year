function [c2,w] = free_neighbours(x0, n)
x = x0;
nodes = [];
w = zeros(1,n);
nodes(1,:) = x0;
w(1)=1;
for i=1:n
    [avail, nmbr] = avail_neigh(x,nodes);
    if nmbr == 0
        break;
    end

    w(i+1)=w(i)*nmbr;
    if nmbr == 1
        xprim = avail;
    else
        xprim = datasample(avail,1);
    end
    nodes(i+1,:) = xprim;
    x=xprim;
end

w=w(2:end);
if length(nodes) - 1 == n
    c2 = 1;
else 
    c2 = 0;
end