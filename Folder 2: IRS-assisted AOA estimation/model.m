%%1. This code is for paper "MUSIC Algorithm for IRS-Assisted AOA Estimation".
%%2. Most of the notations used in the code are consistent with the paper. 
%%3. If you want to obtain the result as in Fig. 4, you need to vary the number of users K.
%%4. If you want to plot the spectrum as in Fig. 3, please use 'plot(search_doa,(P/max(P)),'r')' in the MUSIC function.
%%Author: WANG, Qipeng
%%V1: 2023/08/10
%%V2: 2024/6/13: We update some references for parameter setting, e.g., the
%%antenna gain.
clc
clear all
close all
num_IRS = 128; %%Number of IRS elements
Q = 4;  %%Number of blocks in (18)
K = 3;   %%Number of users
L = 7;  %%Number of elements of the virtual array in (17) and (18)
BW = 10*10^6; %%Bandiwidth of the system
P_t = 10^(23/10)*10^(-3); %Transmit power of each user
MTC = 5000;  %%Number of monte carlo iterations
error = 0;  %%Initialize the number of errors
fc = 2.4; %in GHz
for mtc = 1:MTC
    mtc
    Phi = zeros(L,num_IRS); %%IRS reflection coefficient 
    lambda = 3*10^8/(fc*10^9); %%Wavelength
    d = lambda/2; %%IRS element spacing
    steer_vector_IRS_rfl1 = ((exp(-sqrt(-1)*pi*(0:num_IRS-1)*sin(45*pi/180))).')'; %%Steering vector of the IRS, in this simulation, we assume that the angle between the receiver and the IRS is 45 degree
    B = zeros(num_IRS,K);
    angle_interstore = [];   
 %%%%%%%%Generate user locations
    loc_BS = [0,0];
    loc_IRS = [50,-50];
    d_BI = norm(loc_IRS - loc_BS);
    center_user = [20,-20];
    loc_users = [];
    d_UI = [];
    for k = 1:K
        angle = (unifrnd(0,2))*pi; 
        r = 30*sqrt((rand(1,1)));
        x_k = center_user(1) + r*cos(angle);
        y_k = center_user(2) + r*sin(angle);
        loc_users = [loc_users;x_k,y_k];
        d_UI = [d_UI,norm(loc_users(k,:) - loc_IRS)];
    end
    %%Calculate the AOA from users to the IRS
    source_doa = []; 
    for k = 1:K
        source_doa = [source_doa,-atan((loc_users(k,2)-loc_IRS(2))/(loc_users(k,1)-loc_IRS(1)))*180/pi];
    end
    %%%We need to guarantee that any two users are seperated
    if K >= 2
        ind = nchoosek(1:K,2); 
        angle_diff = abs(source_doa(ind(:,2))-source_doa(ind(:,1))); %difference between any two users
        l = 3;
        while min(angle_diff) < l
            loc_users = [];
            d_UI = [];
            for k = 1:K
                angle = (unifrnd(0,2))*pi; 
                r = 30*sqrt((rand(1,1))); 
                x_k = center_user(1) + r*cos(angle);
                y_k = center_user(2) + r*sin(angle);
                loc_users = [loc_users;x_k,y_k];
                d_UI = [d_UI,norm(loc_users(k,:) - loc_IRS)];
            end
            source_doa = []; 
            for k = 1:K
                source_doa = [source_doa,-atan((loc_users(k,2)-loc_IRS(2))/(loc_users(k,1)-loc_IRS(1)))*180/pi];
            end
            ind = nchoosek(1:K,2); 
            angle_diff = abs(source_doa(ind(:,2))-source_doa(ind(:,1)));
        end 
    end
    %%%calculate the path loss factor, i.e., beta_k in the paper
    beta = [];
    %%%Note that for the transmit antenna gain and the receive antenna gain, we consider the practical antenna product. https://southwestantennas.com/products?filters=omni-antennas&page=2
    G_tx = 6; %antenna gain of the user. The user not necessarily to be a phone, it can be any transmitter, e.g., UAV. According to the above website, the gain of a single antenna can be 6dbi.
    G_rx = 16; %antenna gain of the BS. According to the above website, the antenna gain of the BS can be 16dbi.
    G_irs = 5; %gain of the IRS. Ref: Intelligent Reflecting Surface vs. Decode-and-Forward: How Large Surfaces Are Needed to Beat Relaying?
    for k = 1:K
        beta_k = sqrt(10^((G_tx+G_irs-28-20*log10(fc)-22*log10(d_BI))/10))*sqrt(10^((G_rx+G_irs-28-20*log10(fc)-22*log10(d_UI(k)))/10))*exp(sqrt(-1)*randn(1,1));
        beta = [beta,beta_k]; %%store the path loss factor for each user
    end
    %%%calculate the steering vector at the IRS from K users, i.e., b(theta_k) in (7)
    source_doa = source_doa.*pi/180;  
    for k = 1:K
        angle_interstore = [angle_interstore;source_doa(k)*180/pi];        
        B(:,k) = [(exp(-sqrt(-1)*pi*(0:num_IRS-1)*sin(source_doa(k)))).'];
    end
    angle_store(:,mtc) = angle_interstore;
    %%We use a new notation G_eq to denote \tilde{A} in (21). Note that \tilde{A} = G_eq*B.
    G_eq = [];
    for l = 1:L 
        Phi(l,:) = exp(sqrt(-1)*randn(1,num_IRS));
        G = steer_vector_IRS_rfl1*diag(Phi(l,:));
        G_eq = [G_eq;G];
    end  
    %%Note that X = \bar{A}\tilde{x}is the noise-free signal in (23).
    X = zeros(L,Q);
    virtue_array = [];
    for k =1:K
        virtue_array = [virtue_array,G_eq*B(:,k)];
        X = X + sqrt(P_t)*beta(k)*G_eq*B(:,k)*exp(-sqrt(-1)*2*pi*randn(1,Q)); %%%exp(-sqrt(-1)*2*pi*randn(1, Q)) is the transmitted signal of each user k over Q blocks
    end
    %%Generate the noise
    noise_power = 10^(-169/10)*10^(-3)*BW;
    noise = sqrt(noise_power)*(sqrt(0.5)*randn(L,Q)+sqrt(0.5)*sqrt(-1)*randn(L,Q)); 
    %%Received signal (23)
    Y = X + noise; 
    %%Apply the MUSIC algorithm to estimate angles
    est_angles = MUSIC(Y,Q,K,L,G_eq,num_IRS);
    %%%In the following, we detect the number of errors. If the difference between the true AOA and the estimated AOA is larger than 1 degree, we declare that an error event
    %%%happens.
    true_angles = source_doa*180/pi;
    detect_thre = 1; %%Error detection threshold in degree
    error_mtc = error_detection(est_angles,true_angles,detect_thre,K);
    error = error + error_mtc;
end
%%Calculate the error probability
error_prob = error/MTC/K;
