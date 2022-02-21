function c2 = naive(x0, n)
x = x0;
nodes = [];
nodes(1,:) = x0;
for i=1:n
    pos_dir = posdir_2d(x,1);
    % Draws 1 sample randomly from pos_dir
    xprim = datasample(pos_dir,1);
    % Check if drawn sample has already been visited
    if ismember(xprim,nodes,'rows') == 1
        %ismember(xprim,nodes,'rows')
        break;
    end
    % If not then add new node to visited nodes
    nodes(i+1,:) = xprim;
    x=xprim;
end

if length(nodes) - 1 == n
    c2 = 1;
else 
    c2 = 0;
end
% x1 = -1; y1 = -1;
% x2 = 3;  y2 = 3;
% figure
% hold on
% plot(nodes(:,1),nodes(:,2),'--*')
% axis([x1 x2 y1 y2])
% title("Not a SAW of length n=6")
% xlabel("x")
% ylabel("y")
% grid on
% hold off
% end

