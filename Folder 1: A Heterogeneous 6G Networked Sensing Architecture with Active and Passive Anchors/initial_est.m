function coor_ini=initial_est(distant,coordinate_bs)
%This function aims to initialize targets location estimation. This code is
%taken from "https://github.com/PolyuWirelessCommunicationLab/Device-Free-Sensing-in-OFDM-Cellular-Network/blob/main/initial_est.m".
      a1 = coordinate_bs(1,1);b1=coordinate_bs(1,2);
      a2 = coordinate_bs(2,1);b2=coordinate_bs(2,2);
      a3 = coordinate_bs(3,1);b3=coordinate_bs(3,2);
      r1 = distant(1);
      r2 = distant(2);
      r3 = (distant(3)-r1);
      A = [a1-a2,b1-b2;a1-a3,b1-b3];
      R = [r2^2-r1^2+a1^2+b1^2-a2^2-b2^2;r3^2-r1^2+a1^2+b1^2-a3^2-b3^2]/2;
      coor_ini = A\R;
end