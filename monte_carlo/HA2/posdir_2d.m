function [posx] = posdir_2d(x,n)
%  possible dir for 2d
posx =  zeros(4,2);
posx(1,:) = [x(1)+n,x(2)];
posx(2,:) = [x(1),x(2)+n];
posx(3,:) = [x(1)-n,x(2)];
posx(4,:) = [x(1),x(2)-n];

end