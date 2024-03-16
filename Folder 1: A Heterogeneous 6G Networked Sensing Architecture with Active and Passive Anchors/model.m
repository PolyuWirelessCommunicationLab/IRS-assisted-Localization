clc
clear all
%%%%Parameters zone
%Please vary parameters to obtain the performance under different setups
K_area = [2,2]; %number of targets in each sensing area. It is recommended to make the number of targets be the same in each sensing region.
K = sum(K_area); %number of total targets
x_irs = [-60, 70]; %locations of the IRSs (-60,40) and (70,40)
y_irs = [40, 40]; 
x_bs = [-100, 100]; %locations of the BSs (100,0) and (-100,0)
y_bs = [0, 0];
I = length(x_irs); %%number of IRSs in the network
num_area = I;
l1 = 5;%threshold 1
l2 = 2; %threshold 2
MTC = 10000; %number of monte carlo simulations
BW = 4e8; %bandwidth of the system
c0 = 3e8; %speed of light
%initialize the error number under different detection thresholds
error_num = 0;
for mtc = 1:MTC
    mtc
%%%%randomly generate targets locations and calculate the ranges
    d_D_BS1 = []; %store original direct ranges bewteen targets and BS1, i.e., D_1^{III} in the paper 
    d_D_BS2 = []; %store original direct ranges bewteen targets and BS2, i.e., D_2^{III} in the paper 
    d_R_BS1 = []; %store original indirect ranges bewteen targets and BS1, i.e., D_1^{IV} in the paper  
    d_R_BS2 = []; %store original indirect ranges bewteen targets and BS1, i.e., D_2^{IV} in the paper  
    target_store = [];
    d_BS1_irs = [sqrt((x_bs(1) - x_irs).^2 + (y_bs(1) - y_irs).^2)]; %range between the IRS and BS1
    d_BS2_irs = [sqrt((x_bs(2) - x_irs).^2 + (y_bs(2) - y_irs).^2)]; %range between the IRS and BS2
    d_BS1_BS2 = [sqrt((x_bs(2) - x_bs(1)).^2 + (y_bs(2) - y_bs(1)).^2)];%range between BS1 and BS2
    %%generate targets locations
    for a = 1:num_area
        angle_area = (rand(1,K_area(a)) + 1)*pi;  
        r_area = 50*sqrt((rand(1,K_area(a))));
        x_t_area = x_irs(a) + r_area.*cos(angle_area);
        y_t_area = y_irs(a) - l2 + r_area.*sin(angle_area); %here we need a threshold l1 to avoid the negative range case
        target_store = [target_store;x_t_area',y_t_area']; %store the locations of targets
        %calculate and store all the range information 
        d_D_BS1 = [d_D_BS1,2*sqrt((x_t_area - x_bs(1)).^2 + (y_t_area - y_bs(1)).^2)];
        d_D_BS2 = [d_D_BS2,2*sqrt((x_t_area - x_bs(2)).^2 + (y_t_area - y_bs(2)).^2)];
        d_R_BS1 = [d_R_BS1,d_BS1_irs(a) + sqrt((x_bs(1) - x_t_area).^2 + (y_bs(1) - y_t_area).^2) + sqrt((x_irs(a) - x_t_area).^2 + (y_irs(a) - y_t_area).^2)];
        d_R_BS2 = [d_R_BS2,d_BS2_irs(a) + sqrt((x_bs(2) - x_t_area).^2 + (y_bs(2) - y_t_area).^2) + sqrt((x_irs(a) - x_t_area).^2 + (y_irs(a) - y_t_area).^2)]; 
    end
    %calculate the delay information in terms of OFDM samples
    tau_D_BS1 = ceil(d_D_BS1*BW/c0);
    tau_D_BS2 = ceil(d_D_BS2*BW/c0);
    tau_R_BS1 = ceil(d_R_BS1*BW/c0);
    tau_R_BS2 = ceil(d_R_BS2*BW/c0);
    %calculate the estimated range information using OFDM samples, i.e., \hat{D}'s in the paper  
    d_D_BS1 = ((tau_D_BS1 - 0.5)*c0/2/BW);
    d_D_BS2 = ((tau_D_BS2 - 0.5)*c0/2/BW);
    d_R_BS1 = (tau_R_BS1 - 0.5)*c0/BW;
    d_R_BS2 = (tau_R_BS2 - 0.5)*c0/BW; 
    %Note that we need to guarantee that different targets are in different OFDM delay samples. If this condition does not hold, we generate new targets locations
    while (length(unique([tau_D_BS1])) ~= K)||(length(unique([tau_D_BS2])) ~= K)||(length(unique([tau_R_BS1])) ~= K)||(length(unique([tau_R_BS2])) ~= K)||(min(abs(d_D_BS1 + d_D_BS2 - d_BS1_BS2)) < l1)||(min(abs(max(d_D_BS1,d_D_BS2) -min(d_D_BS1,d_D_BS2) - d_BS1_BS2)) < l1)
        d_D_BS1 = [];
        d_D_BS2 = [];
        d_R_BS1 = [];
        d_R_BS2 = [];
        target_store = [];
        d_BS1_irs = [sqrt((x_bs(1) - x_irs).^2 + (y_bs(1) - y_irs).^2)];
        d_BS2_irs = [sqrt((x_bs(2) - x_irs).^2 + (y_bs(2) - y_irs).^2)];
        d_BS1_BS2 = [sqrt((x_bs(2) - x_bs(1)).^2 + (y_bs(2) - y_bs(1)).^2)];
        for a = 1:num_area
            angle_area = (rand(1,K_area(a))+1)*pi; 
            r_area = 50*sqrt((rand(1,K_area(a))));
            x_t_area = x_irs(a) + r_area.*cos(angle_area);
            y_t_area = y_irs(a) - l1 + r_area.*sin(angle_area);
            target_store = [target_store;x_t_area',y_t_area'];
            d_D_BS1 = [d_D_BS1,sqrt((x_bs(1) - x_t_area).^2 + (y_bs(1) - y_t_area).^2)];
            d_D_BS2 = [d_D_BS2,sqrt((x_bs(2) - x_t_area).^2 + (y_bs(2) - y_t_area).^2)];
            d_R_BS1 = [d_R_BS1,d_BS1_irs(a) + sqrt((x_bs(1) - x_t_area).^2 + (y_bs(1) - y_t_area).^2) + sqrt((x_irs(a) - x_t_area).^2 + (y_irs(a) - y_t_area).^2)];
            d_R_BS2 = [d_R_BS2,d_BS2_irs(a) + sqrt((x_bs(2) - x_t_area).^2 + (y_bs(2) - y_t_area).^2) + sqrt((x_irs(a) - x_t_area).^2 + (y_irs(a) - y_t_area).^2)]; 
        end
        tau_D_BS1 = ceil(2*d_D_BS1*BW/c0);
        tau_D_BS2 = ceil(2*d_D_BS2*BW/c0);
        tau_R_BS1 = ceil(d_R_BS1*BW/c0);
        tau_R_BS2 = ceil(d_R_BS2*BW/c0);      
        %calculate the estimated range information using OFDM samples, i.e., \hat{D}'s in the paper  
        d_D_BS1 = ((tau_D_BS1 - 0.5)*c0/2/BW);
        d_D_BS2 = ((tau_D_BS2 - 0.5)*c0/2/BW);
        d_R_BS1 = (tau_R_BS1 - 0.5)*c0/BW;
        d_R_BS2 = (tau_R_BS2 - 0.5)*c0/BW;    
    end

    tic
    %%%%In the following, we perform data association and localization based on \hat{D}'s in the above.
    %First, we find the data association for \gamma_k's based on problem (40)
    da_phaseI = [];
    for i = 1:K
        for j = 1:K
            if d_D_BS1(i) + d_D_BS2(j) > d_BS1_BS2 - l2 && abs(d_D_BS1(i)-d_D_BS2(j)) < d_BS1_BS2 + l2 
                if I >= 2 %In the multi-IRS case, we need to associated \gamma_k's
                    [x_pssbl,y_pssbl] = sol_two_equa(x_bs,y_bs,d_D_BS1(i),d_D_BS2(j));           
                    if length(x_pssbl) ~= 0 && sum(isnan(x_pssbl)) == 0
                        x_pssbl = double(x_pssbl);
                        y_pssbl = double(y_pssbl);
                        for x_p = 1:2
                            d_I = zeros(1,num_area);
                                for a = 1:num_area
                                    d_I(a) = sqrt((x_irs(a) - x_pssbl(x_p))^2 + (y_irs(a) - y_pssbl(x_p))^2);
                                end
                            da_phaseI = [da_phaseI;i,j,find(d_I == min(d_I)),min(d_I)];
                        end
                    end
                end
                if I == 1 %In the single IRS case, we don't need to do that.
                   da_phaseI = [da_phaseI;i,j,1,0];
                end
            end    
        end
    end
    da_phaseI = unique(da_phaseI,'rows','stable');
    [e_row,e_col] = size(da_phaseI);
    da_phaseII = [];
    da_phaseII_B1 = [];
    da_phaseII_B2 = [];
    %Second, according to (24) and (25), all the ranges should be positive.
    for k1 = 1:K
        for e = 1:e_row
            if I >= 2 %multi-IRS case
                if d_R_BS1(k1) - d_D_BS1(da_phaseI(e,1)) - d_BS1_irs(da_phaseI(e,3)) > 0
                    da_phaseII_B1 = [da_phaseII_B1;k1,da_phaseI(e,1),da_phaseI(e,3),d_R_BS1(k1) - d_D_BS1(da_phaseI(e,1)) - d_BS1_irs(da_phaseI(e,3))];
                end
            end
            if I == 1 %single_IRS case
                if d_R_BS1(k1) - d_D_BS1(da_phaseI(e,1)) - d_BS1_irs > 0
                    da_phaseII_B1 = [da_phaseII_B1;k1,da_phaseI(e,1),da_phaseI(e,3),d_R_BS1(k1) - d_D_BS1(da_phaseI(e,1)) - d_BS1_irs];
                end
            end
        end
    end
    da_phaseII_B1 = unique(da_phaseII_B1,'rows','stable');
    for k2 = 1:K
        for e2 = 1:e_row
            if I >= 2
                if d_R_BS2(k2) - d_D_BS2(da_phaseI(e2,1)) - d_BS2_irs(da_phaseI(e2,3)) > 0
                    da_phaseII_B2 = [da_phaseII_B2;k2,da_phaseI(e2,1),da_phaseI(e2,3),d_R_BS2(k2) - d_D_BS2(da_phaseI(e2,2)) - d_BS2_irs(da_phaseI(e2,3))];
                end
            end
            if I == 1
                if d_R_BS2(k2) - d_D_BS2(da_phaseI(e2,1)) - d_BS2_irs > 0
                    da_phaseII_B2 = [da_phaseII_B2;k2,da_phaseI(e2,1),da_phaseI(e2,3),d_R_BS2(k2) - d_D_BS2(da_phaseI(e2,2)) - d_BS2_irs];
                end
            end
        end
    end
    da_phaseII_B2 = unique(da_phaseII_B2,'rows','stable');
    [da_II_B1_row,da_II_B1_col] = size(da_phaseII_B1);
    [da_II_B2_row,da_II_B2_col] = size(da_phaseII_B2);
    %Third, we find data associations that satisfy (34).
    for b1 = 1:da_II_B1_row
        for b2 = 1:da_II_B2_row
            if abs(da_phaseII_B1(b1,4) - da_phaseII_B2(b2,4)) < l2
                da_II = [da_phaseII_B1(b1,1:2),da_phaseII_B2(b2,1:3)];
                num = zeros(1,K);
                for k =1:K
                    num(1,k) = sum(da_II(1:4) == k);
                end
                da_phaseII = [da_phaseII;da_II,num];
            end   
        end
    end
    da_phaseII = unique(da_phaseII,'rows','stable');
    [da_II_row,da_II_col] = size(da_phaseII);
    for k = 1:K
        for i = 1:da_II_row
            if da_phaseII(i,1) == k
                r(k+1) = i;
            end
        end
    end
    %%Repeat iterations in Algorithm 1 to remove data assocaitions that does not satisfy (37)
    est_diff_store = [];
    %%estimate targets locations given any data association solution
    for d = 1:r(K+1)
          dist = [d_D_BS1(da_phaseII(d,2)),d_D_BS2(da_phaseII(d,4)),d_R_BS1(da_phaseII(d,1)) - d_BS1_irs(da_phaseII(d,5)),d_R_BS2(da_phaseII(d,3)) - d_BS2_irs(da_phaseII(d,5))];
          true_dis = [d_D_BS1(da_phaseII(d,2)),d_D_BS2(da_phaseII(d,4)),d_R_BS1(da_phaseII(d,1)),d_R_BS2(da_phaseII(d,3))];
          coordinate_bs = [x_bs(1),y_bs(1);
                        x_bs(2),y_bs(2);
                        x_irs(da_phaseII(d,5)),y_irs(da_phaseII(d,5));
                        x_irs(da_phaseII(d,5)),y_irs(da_phaseII(d,5))]; 
           [diff2,coordinate_estimation,diff1] = localization_with_DA(dist,coordinate_bs); 
           x_est = coordinate_estimation(1);
           y_est = coordinate_estimation(2);  
           est_d_D_BS1 = sqrt((x_est - x_bs(1))^2+(y_est - y_bs(1))^2);
           est_d_D_BS2 = sqrt((x_est - x_bs(2))^2+(y_est - y_bs(2))^2);
           est_d_R_BS1 = sqrt((x_est - x_bs(1))^2+(y_est - y_bs(1))^2) + sqrt((x_est - x_irs(da_phaseII(d,5)))^2 + (y_est - y_irs(da_phaseII(d,5)))^2)+ sqrt((x_bs(1) - x_irs(da_phaseII(d,5)))^2 + (y_bs(1) - y_irs(da_phaseII(d,5)))^2);
           est_d_R_BS2 = sqrt((x_est - x_bs(2))^2+(y_est - y_bs(2))^2) + sqrt((x_est - x_irs(da_phaseII(d,5)))^2 + (y_est - y_irs(da_phaseII(d,5)))^2)+ sqrt((x_bs(2) - x_irs(da_phaseII(d,5)))^2 + (y_bs(2) - y_irs(da_phaseII(d,5)))^2);         
           est_dis = [est_d_D_BS1,est_d_D_BS2,est_d_R_BS1,est_d_R_BS2];
           est_diff_store = [est_diff_store;norm(est_dis - true_dis)^2];
           cost(d) = norm(est_dis - true_dis)^2; 
        end     
        da_crrct = [];
        da_crrct_store = [];
        for d = 1:r(K+1)
            if cost(d) <= 2
               da_crrct = [da_crrct;da_phaseII(d,1:5),cost(d)];
            end
        end
        [da_crrct_row,da_crrct_col] = size(da_crrct);
        %%Calculate the objective values of different data association solution and find the optimal one
        if da_crrct_row ~= K
            da_final = checkwithconstraints(K,da_crrct,da_crrct_row);
            da_crrct = da_final;
        end
        [da_crrct_row,da_crrct_col] = size(da_crrct);
        cost_overall = [];
        if da_crrct_row ~= K
            overall_group = da_crrct_row/K;
            for g = 1:overall_group
                cost_overall(g) = sum(da_final([K*(g-1)+1:K*g],6));
            end
            ind_min = find(cost_overall == min(cost_overall));
            da_crrct = da_final([K*(ind_min-1)+1:K*ind_min],:);
        end
        time(mtc) = toc;
        if ~isempty(da_crrct)
            est_target = zeros(K,2);
            %%estimate targets locations given any data association solution
            for d = 1:K
                  dist = [d_D_BS1(da_crrct(d,2)),d_D_BS2(da_crrct(d,4)),d_R_BS1(da_crrct(d,1)) - d_BS1_irs(da_crrct(d,5)),d_R_BS2(da_crrct(d,3)) - d_BS2_irs(da_crrct(d,5))];
                  true_dis = [d_D_BS1(da_crrct(d,2)),d_D_BS2(da_crrct(d,4)),d_R_BS1(da_crrct(d,1)),d_R_BS2(da_crrct(d,3))];
                  coordinate_bs = [x_bs(1),y_bs(1);
                                    x_bs(2),y_bs(2);
                                    x_irs(da_crrct(d,5)),y_irs(da_crrct(d,5));
                                    x_irs(da_crrct(d,5)),y_irs(da_crrct(d,5))]; 
                   [diff,coordinate_estimation,diff1] = localization_with_DA(dist,coordinate_bs); 
                   x_est = coordinate_estimation(1);
                   y_est = coordinate_estimation(2);
                   est_target(d,:) = [x_est,y_est];
            end
            %%perform error detection under different detection thresholds
            [diff,rowindex] = optimal_assign(est_target,target_store,K);
            error_num = error_num + length(find(diff > 0.8));
        end
end
%calculate the error probability
error_prob = error_num/K/MTC;


