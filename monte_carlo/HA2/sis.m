function [X, omegas] = sis(X_old,omega_old)
    [r,c,z] = size(X_old);
    X = zeros(r+1,c,z);
    X(1:r,:,:) = X_old;
    omegas = ones(z,2);
    omegas(:,1) = omega_old(:,1);
    if r == 1
         Xn = X(1,:,:);
    else
         Xn = X(end-1,:,:);
    end

    for i=1:z
        if r == 1
            [avail,nbr] = avail_neigh(Xn(1,:,i), X(1,:,i));
        else
            [avail,nbr] = avail_neigh(Xn(1,:,i), X(1:end-1,:,i));
        end
        
        if nbr == 0
            xprim = zeros(c,1);
        elseif nbr == 1
            xprim = avail;
        else
            xprim = datasample(avail,1);
        end
        X(end,:,i) = xprim;
        omegas(i,2) = omegas(i,1)*nbr;
    end
    %omega_old
    %omegas
end