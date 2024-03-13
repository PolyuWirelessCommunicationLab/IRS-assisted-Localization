function error_mtc = error_detection(est_angles,true_angles,detect_thre,K)
%error_detection detects errors based on the estimated angles and true
%angles. The error detection problem is formulated as a Linear Assignment
%Problem and solved using the 'munkres' function.

%Inputs: 
%   est_angles: the estimated angles
%   true_angles: the true angles
%   detect_thre: detection thereshold
%   K: the number of users
%Output: 
%   error_mtc: the number of detected errors
    diff = zeros(1,K);
    diffdistant = zeros(K,K);
    for i = 1:K
        for j = 1:K
            diffdistant(i,j) = abs(est_angles(i) - true_angles(j));
        end
    end
    [rowindex,cost] = munkres(diffdistant);
    cor_index = rowindex;
    for i = 1:K
        diff(i) = diffdistant(i,rowindex(i));  
    end
    error_mtc = length(find(diff>detect_thre));
end