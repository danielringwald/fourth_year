function [posx] = posdir_multi(x,d)
%Gives the availablie direction from x with stepsize 1 with dimensions d
posx = repelem(x,2*length(x),1);

for ii=1:d
    posx(ii,ii) = x(ii) + 1;
    posx(ii+d,ii) = x(ii) - 1;
end
end
