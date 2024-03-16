function [diff,coordinate_estimation,diff1] = localization_with_DA(dist,coordinate_bs,coordinate_tg)
    %This function aims to solve problem (30) in an iterative manner. This
    %code is modified from "https://github.com/PolyuWirelessCommunicationLab/Device-Free-Sensing-in-OFDM-Cellular-Network/blob/main/coorest.m".
    l_d = length(dist);
    coorini = initial_est(dist,coordinate_bs); 
    x_est = coorini(1);
    y_est = coorini(2);
    dist(3) = (dist(3)-dist(1));
    dist(4) = (dist(4)-dist(2));
    beta_current = [x_est;y_est];
    beta_next = zeros(2,1);
    J = zeros(l_d,2);
    sdiff = 10;
    itera = 1;
    sdiff_next = 0;
    sdiff_current = 0;
    while sdiff>0.02 && itera < 100
        x_est=beta_current(1);
        y_est=beta_current(2);
        r=zeros(l_d,1);
        for i = 1:4
            J(i,1) = (x_est-coordinate_bs(i,1))/sqrt((x_est-coordinate_bs(i,1))^2+(y_est-coordinate_bs(i,2))^2);
            J(i,2) = (y_est-coordinate_bs(i,2))/sqrt((x_est-coordinate_bs(i,1))^2+(y_est-coordinate_bs(i,2))^2);
            r(i,1) = sqrt((x_est-coordinate_bs(i,1))^2+(y_est-coordinate_bs(i,2))^2)-dist(i);
        end 
        lamada = 0.2;
        beta_next = beta_current - lamada*pinv(J.'*J)*J.'*r;
        sdiff = norm(beta_next - beta_current);
        itera = itera + 1;
        beta_current = beta_next;
    end
    coordinate_estimation = beta_current.';
    diff = [];
    for n = 1:l_d
        diff(n) = abs(norm(coordinate_estimation - coordinate_bs(n,:)) - dist(n));
    end
    diff = sum(diff);
    diff1 = [];
    for n = 1:l_d
        diff1(n) = abs(norm(coorini' - coordinate_bs(n,:)) - dist(n));
    end
    diff1 = (sum(diff1));
end