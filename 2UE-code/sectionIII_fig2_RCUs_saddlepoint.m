function sectionIII_fig2_RCUs_saddlepoint(n, np, rho_db, b, angle2, Mlist, nbrOfRealizations, ESTIMATOR, COMBINER, PILOT_CONTAMINATION, UNCORRELATED, NO_OF_UEs, BS_CSI, UE_CSI)
% Function sectionIII_fig2_RCUs_saddlepoint(n, np, rho_db, b, angle2, Mlist, nbrOfRealizations, ESTIMATOR, COMBINER, PILOT_CONTAMINATION, UNCORRELATED, NO_OF_UEs, BS_CSI, UE_CSI): 
% Generates saddlepoint approximations in Fig. 2 for both UL and DL. 
% 
% INPUTS:
% n = blocklength
% np = length pilot sequence
% rho_db = transmit power [dB]
% b = information bits
% angle_2 = angle of UE2 to BS
% Mlist = list with number of antennas at BS
% nbrOfRealizations = number of realization in the Monte-Carlo simulations
% ESTIMATOR = Estimator to use [LS, MMSE]
% COMBINER = Combiner to use [MR, M-MMSE, RZF]
% PILOT_CONTAMINATION = [0,1]; 1: the two UEs use the same pilot sequence
% UNCORRELATED = [0,1]; 0: spatially correlated channel
% NO_OF_UEs = number of UEs [1,2]
% BS_CSI = [0,1]; 1: Perfect CSI at BS
% UE_CSI = [0,1]; 1: Perfect CSI at UE
%
 
DEBUG = 1;

if DEBUG == 1
    %FLAGS
    np = 2; %number of pilots
    n = 300; %total number of channle uses
    b = 8*20; %information bits
    rho_db = 10; %transmit power [dB]
    angle2 = deg2rad(0:5:65); % angle of UE2 to BS
    nbrOfRealizations=1e5; %number of saddlepoint realizations
    Mlist = [100]; %the number of antennas considered
    COMBINER = 'MR'; %what combiner to use [MR, M-MMSE, RZF]
    ESTIMATOR = 'MMSE'; %what estimator to use [LS, MMSE]
    PILOT_CONTAMINATION = 1; %Let the two UEs use the same pilot sequence
    NO_OF_UEs = 2; % 1 or 2
    UNCORRELATED = 0; %0 = spatial correlation
    BS_CSI = 0; %1 => Perfect CSI at BS
    UE_CSI = 0; %0 = UE uses E[H1'W1]; 1 = UE knows H1'W1; 2= UE knows W1 and Hhat1
end
% Initializations:
angle1 = deg2rad(30); % angle of UE1 to BS
ASDdeg = 25; %angular delay spread

n_ul = (n-np)/2; % available channel uses for the UL
n_dl = n_ul; % available channel uses for the DL
rho = 10^(rho_db/10);  %transmit power [linear]

format long
avg_error_ul = inf(length(angle2),length(Mlist));
avg_error_dl = inf(length(angle2),length(Mlist));
avg_error = inf(length(angle2),length(Mlist), 2);
s_val_ul = nan(length(angle2), length(Mlist));
s_val_dl=nan(length(angle2), length(Mlist));

% For each angle:
for j = 1:length(angle2)
    
    M=Mlist; %consider current M
    
    [H1, H2, R1, R2] = generateChannel(M, angle1, angle2(j), ASDdeg, nbrOfRealizations,UNCORRELATED);
    
    if BS_CSI == 0
        [hhat1, hhat2, C1, C2] =  estimateChannel(rho, np, H1, H2, R1, R2, PILOT_CONTAMINATION, ESTIMATOR);      
    else
        hhat1 = H1; C1 = zeros(M,M); hhat2 = H2; C2 = zeros(M,M); %Perfect CSI at BS
    end
    
    disp('Combiner')
    tic
    [v1, v2] = createCombiners(rho, hhat1, hhat2, C1, C2, COMBINER, NO_OF_UEs);
    toc
    %-----------------------------------------------
    % UPLINK
    %-----------------------------------------------
    % Initializations:
    sigma_sq_list = nan(1,nbrOfRealizations);
    g_list = nan(1,nbrOfRealizations);
    g_hat_list= nan(1,nbrOfRealizations);
    % Computation of effective noise, effective channel and channel est.
    for indx = 1:nbrOfRealizations
        v1_cur =  v1(:,indx);
        if NO_OF_UEs == 1
            sigma_sq_list(indx) = v1_cur'*v1_cur;
        elseif NO_OF_UEs==2
            sigma_sq_list(indx) = rho*abs(v1_cur' * H2(:,indx))^2 + v1_cur'*v1_cur;
        end
        g_list(indx) = v1_cur' * H1(:,indx) ;
        g_hat_list(indx) = v1_cur' * hhat1(:,indx);
    end
    % Computation of the UL average error probability using the saddlepoint approx:
    [avg_error_ul(j), s_val_ul(j)] = sec3_fig2_getErrorProbability(n_ul, rho, b, g_list, g_hat_list, sigma_sq_list, nbrOfRealizations, 1);
    disp(['UL error probability for (angle2, M) = (' num2str(rad2deg(angle2(j))) ',' num2str(M) '): ' num2str(avg_error_ul(j)) ' and s = ' num2str(s_val_ul(j))]);
    toc
    %-----------------------------------------------
    % DOWNLINK
    %-----------------------------------------------
    tic
    w1 = v1 ; %precoding vector for UE1
    w2 = v2 ; %precoding vector for UE2
    % Initializations:
    sigma_sq_list_dl = nan(1,nbrOfRealizations);
    g_list_dl = nan(1,nbrOfRealizations);
    ghat_dl_list = nan(1,nbrOfRealizations);
    % Computation of effective channel and channel estimation:
    if UE_CSI == 0
        ghat_dl = mean(sum(conj(H1).*w1,1)); %UE uses channel hardening
        ghat_dl_list = ghat_dl*ones(1,nbrOfRealizations);
    elseif UE_CSI == 1
        ghat_dl_list = sum(conj(H1).*w1,1); %UE knows the precoded channel
    elseif UE_CSI == 2
        ghat_dl_list = sum(conj(hhat1).*w1,1); %UE knows Hhat1 and assumes it to be correct
    end
    % Computation of effective noise:
    for indx = 1:nbrOfRealizations
        if NO_OF_UEs == 1
            sigma_sq_list_dl(indx) =  1;
        elseif NO_OF_UEs==2
            sigma_sq_list_dl(indx) = rho*abs(H1(:,indx)'*w2(:,indx))^2 + 1;
        end
        g_list_dl(indx) = H1(:,indx)'*w1(:,indx);
    end
    % Computation of the DL average error probability:
    [avg_error_dl(j) ,  s_val_dl(j)] = sec3_fig2_getErrorProbability(n_dl, rho, b, g_list_dl, ghat_dl_list, sigma_sq_list_dl, nbrOfRealizations, 1);
    disp(['DL error probability for (angle2, M) = (' num2str(rad2deg(angle2(j))) ',' num2str(M) '): ' num2str(avg_error_dl(j)) ' and s = ' num2str(s_val_dl(j))]);
    toc
end
if DEBUG == 1
    
    figure(1); hold on;
    semilogy(rad2deg(angle2), avg_error_ul) ;
    set(gca, 'YScale', 'log')
    xlabel('Angle of UE 2')
    title('UL')
    ylim([1e-5,1])
    
    figure(2); hold on;
    semilogy(rad2deg(angle2), avg_error_dl) ;
    set(gca, 'YScale', 'log')
    xlabel('Angle of UE 2')
    title('DL')
    ylim([1e-5,1])
    
end
avg_error(:,:,1) = avg_error_ul;
avg_error(:,:,2) = avg_error_dl;

avg_error(isinf(avg_error)) = nan;

%--------------------------------------------------------------
%   Save file
%--------------------------------------------------------------
data.np=np;
data.n = n;
data.rate =   b / n_ul;
data.snr_db = rho_db;
data.angle1 = angle1;
data.angle2 = angle2;
data.nbrOfRealizations = nbrOfRealizations;
data.Mlist = Mlist;
data.COMBINER = COMBINER;
data.ESTIMATOR = ESTIMATOR;
data.PILOT_CONTAMINATION=PILOT_CONTAMINATION;
data.NO_OF_UEs = NO_OF_UEs;
data.UNCORRELATED = UNCORRELATED;
data.s_val_ul = s_val_ul;
data.avg_error_ul = avg_error(:,:,1);
data.s_val_dl = s_val_dl;
data.avg_error_dl = avg_error(:,:,2);

filename = [ESTIMATOR '_' COMBINER '_UEs_' num2str(NO_OF_UEs) '_SNR_' num2str(rho_db) '_np_' num2str(np) '_n_' num2str(n) '_M_' num2str(Mlist) '.mat'];


if UNCORRELATED == 1
    filename = ['Uncorrelated_' filename];
else
    filename = ['SpatialCorrelation_' filename];
end

if PILOT_CONTAMINATION == 1
    filename = ['pilot_contamination_' filename];
end

if BS_CSI == 1
    filename = ['BS_CSI_' filename];
end

if UE_CSI == 1
    filename = ['UE_CSI_' filename];
elseif UE_CSI == 2
    filename = ['UE_iCSI_' filename];
end
filename = ['RCUs_SP_ULA_Sweep_Angle_' filename];
save(filename, 'data', '-v7.3');


end


function [H1, H2, R1, R2] = generateChannel(M, angle1, angle2, ASDdeg, nbrOfRealizations, uncorr)
% Function to generate the channels according to the book "Massive MIMO
% Networks" by E. Björnson, J. Hoydis and L. Sanguinetti. 

%---------------------------------
% Obtain covariance matrices
%---------------------------------
if uncorr == 0
    R1 = functionRlocalscattering(M,angle1,ASDdeg,1/2); %Covariance matrix of UE 1
    R2 = functionRlocalscattering(M,angle2,ASDdeg,1/2); %Covariance matrix of UE 2
    
    %Circular array
    %     R1 = functionRlocalscatteringApproxUCA(M,angle1,ASDdeg);
    %     R2 = functionRlocalscatteringApproxUCA(M,angle2,ASDdeg);
else
    R1 = eye(M);
    R2 = eye(M);
end
%---------------------------------
%Create channel matrices
%---------------------------------
H1 = sqrt(0.5)*sqrtm(R1)*(randn(M,nbrOfRealizations)+1i*randn(M,nbrOfRealizations));
H2 = sqrt(0.5)*sqrtm(R2)*(randn(M,nbrOfRealizations)+1i*randn(M,nbrOfRealizations));

end

function [hhat1, hhat2, C1, C2] =  estimateChannel(rho, np, H1, H2, R1, R2, PILOT_CONTAMINATION, ESTIMATOR)
% Function to estimate the channels according to the book "Massive MIMO
% Networks" by E. Björnson, J. Hoydis and L. Sanguinetti. 

sigma_ul = 1;
M = size(H1,1);
nbrOfRealizations = size(H1,2);
%---------------------------------
% Estimate the channel for UE 1 and UE 2 to the BS serving UE 1
%---------------------------------

Np = sqrt(0.5)*(randn(M,nbrOfRealizations) + 1i*randn(M,nbrOfRealizations)); %noise on the channel estimation

if PILOT_CONTAMINATION == 1
    yp1 = rho*np*H1 + rho*np*H2 + sqrt(np*rho)*Np; %UE1 and UE 2 use the same pilot sequence
    yp2=yp1;
    Q1 = (rho*np*R1 + rho*np*R2 + sigma_ul*eye(M)); % matrix for MMSE estimation (common for UE1 and UE2)
    Q2=Q1;
else
    yp1 = rho*np*H1 + sqrt(np*rho)*Np;  % UE 1 and UE2 does not use the same pilot sequence
    yp2 = rho*np*H2 + sqrt(np*rho)*Np;  % UE 1 and UE2 does not use the same pilot sequence
    Q1 = (rho*np*R1 + sigma_ul*eye(M)); % matrix for MMSE estimation of H1
    Q2 = (rho*np*R2 + sigma_ul*eye(M)); % matrix for MMSE estimation of H2
end

if strcmp(ESTIMATOR, 'MMSE')
    R1Qinv = R1 / Q1;
    R2Qinv = R2 / Q2;
    
    C1 = R1 - np*rho*R1Qinv*R1; %covariance of Hhat1
    C2 = R2 - np*rho*R2Qinv*R2; %covariance of Hhat2
    
    hhat1 = R1Qinv*yp1; %estimate H1
    hhat2 = R2Qinv*yp2; %estimate H2
elseif strcmp(ESTIMATOR, 'LS')
    A_LS = 1/(rho*np);
    hhat1 = A_LS*yp1; %estimate H1
    hhat2 = A_LS*yp2; %estimate H2    
    % The following can be wrongly normalized, but it is not used anyways
    productAR1 = sqrt(rho)*A_LS * R1; % This sqrt(rho) is introduced to normalize as in Luca's code. Maybe wrong.
    productAR2 = sqrt(rho)*A_LS * R2;
    C1 = R1 - (productAR1+productAR1')*sqrt(rho)*np+np*A_LS*Q1*A_LS'; %covariance of Hhat1
    C2 = R2 - (productAR2+productAR2')*sqrt(rho)*np+np*A_LS*Q2*A_LS'; %covariance of Hhat2
end
end
function [v1, v2] = createCombiners(rho, hhat1, hhat2, C1, C2, COMBINER, NO_OF_UEs)
% Function to create the combiners/precoders according to the book "Massive MIMO
% Networks" by E. Björnson, J. Hoydis and L. Sanguinetti. 

vecnorm = @(A)  sqrt(sum(A.*conj(A),1)); %norm of each column in matrix
sigma_ul = 1;
M = size(hhat1,1);
nbrOfRealizations = size(hhat1,2);
%---------------------------------
%   Create Combining vector
%---------------------------------
if strcmp(COMBINER, 'MR')
    v1 = hhat1 ./ (vecnorm(hhat1) ); %MR combiner for UE 1
    v2 = hhat2 ./ (vecnorm(hhat2) ); %MR combiner for UE 1

elseif strcmp(COMBINER, 'RZF') 
    V_RZF = nan(M,nbrOfRealizations,NO_OF_UEs);
    for k = 1:nbrOfRealizations %estimate each of the realizations
        V = [hhat1(:,k) hhat2(:,k)];
        V_RZF(:,k,:) = rho*V/(rho*(V'*V)+eye(2)); %RZF combiner
    end
    v1 = V_RZF(:,:,1);
    v2 = V_RZF(:,:,2);
    v1 = v1./vecnorm(v1);
    v2 = v2./vecnorm(v2);
    
elseif strcmp(COMBINER, 'M-MMSE')
    
    v1 = nan(size(hhat1));
    v2 = nan(size(hhat2));
    if NO_OF_UEs == 1
        Z = C1 + sigma_ul*eye(M); %matrix used for the inverse in M-MMSE estimation
        
        parfor k = 1:nbrOfRealizations %estimate each of the realizations
            v1(:,k) = rho* ( (rho*(hhat1(:,k)*hhat1(:,k)') + Z) \ hhat1(:,k)); % M-MMSE combiner for UE 1
            v1(:,k) = v1(:,k)/norm(v1(:,k));
        end
        
    elseif NO_OF_UEs == 2
        Z = C1 + C2 + sigma_ul*eye(M); %matrix used for the inverse in M-MMSE estimation
        
        parfor k = 1:nbrOfRealizations %estimate each of the realizations
            PHI = (rho*(hhat1(:,k)*hhat1(:,k)') + rho*(hhat2(:,k)*hhat2(:,k)') + Z);
            
            v1(:,k) = rho* ( PHI \ hhat1(:,k)); % M-MMSE combiner for UE 1
            v1(:,k) = v1(:,k)/norm(v1(:,k));
            
            v2(:,k) = rho * ( PHI \ hhat2(:,k)); % M-MMSE combiner for UE 2
            v2(:,k) = v2(:,k)/norm(v2(:,k));
        end
        
    end
end
end
