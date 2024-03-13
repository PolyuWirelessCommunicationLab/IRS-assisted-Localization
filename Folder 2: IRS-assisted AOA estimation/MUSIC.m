function est_angles = MUSIC(Y,Q,K,L,G_eq,num_IRS)
%MUSIC Estimates AOAs based on the created multi-dimensional signals.
%Inputs: 
%   Y: the received signal
%   Q: the number of blocks
%   K: the number of users
%   L: the number of virtual array elements
%   num_IRS: the number of IRS elements
%Output: 
%   est_angles: The estimated angles
    R = Y*Y';
    R = R/Q;
    [EVector,EValue] = eig(R);
    EVA = diag(EValue);
    [EVA,I] = sort(EVA);
    EVA = fliplr(EVA);
    EVector = fliplr(EVector(:,I));
    search_doa = 0:0.05:90;  %%Feel free to change this
    P = []; %%Store the
    for i=1:length(search_doa)
        a = G_eq*exp(-sqrt(-1)*pi*(0:num_IRS-1)*sin(search_doa(i)*pi/180)).';
        a = a'/norm(a);
        EN = EVector(:,(K+1):L);
        P = [P, abs(1./(a*EN*EN'*a'))];
    end
    %plot(search_doa,(P/max(P)),'r'); %plot the spectrum, e.g., Fig. 3 in the paper  
    [est_angles, est_idx, resolved] = find_doa_from_spectrum_1d(search_doa, P, K); %find peaks
end