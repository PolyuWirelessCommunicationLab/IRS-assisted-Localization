function [diff,rowindex] = optimal_assign(A,B,K)
%This function aims to associate the estimated target location with the true target location.
    diff = zeros(1,K);
    diffdistant = zeros(K,K);
    for i = 1:K
        for j = 1:K
            ab = A(i,:) - B(j,:);
            diffdistant(i,j) = norm(ab);
        end
    end
    [rowindex,cost] = munkres(diffdistant);
    cor_index = rowindex;
    for i = 1:K
      diff(i) = diffdistant(i,rowindex(i));  
    end
end