function sectionIII_fig4_more_options(n, np, rho_db, b, angle2, Mlist, nbrOfRealizations, ESTIMATOR, COMBINER, PILOT_CONTAMINATION, UNCORRELATED, NO_OF_UEs, BS_CSI, UE_CSI, antennaType)
% Function sectionIII_fig4(n, np, rho_db, rate, angle2, Mlist, nbrOfRealizations, COMBINER, PILOT_CONTAMINATION, PERFECT_CSI, NO_OF_UEs, OUTAGE): 
% Generates saddlepoint approximations in Fig. 2 for both UL and DL. 
% 
% INPUTS:
% n = blocklength
% np = length pilot sequence
% rho_db = transmit power [dB]
% b = information bits
% angle2 = angle UE2 to BS
% Mlist = list with number of antennas at BS
% nbrOfRealizations = number of realization in the Monte-Carlo simulations
% ESTIMATOR = Estimator to use [LS, MMSE]
% COMBINER = Combiner to use [MR, M-MMSE, RZF]
% PILOT_CONTAMINATION = [0,1]; 1: the two UEs use the same pilot sequence
% UNCORRELATED = [0,1]; 0: spatially correlated channel
% NO_OF_UEs = number of UEs [1,2]
% BS_CSI = [0,1]; 1: Perfect CSI at BS
% UE_CSI = [0,1,2]; 0: channel hardening; 1: Perfect CSI at UE; 2: UE knows
%                   precoder and channel estimates (same info as in the UL at the BS)
% antennaType = ['ULA','UCA'] uniform linear array or uniform circular array
DEBUG = 1;

if DEBUG == 1
    %FLAGS
    np = 2; %number of pilots
    n = 300; %total number of channle uses
    b = 8*20; %transmission rate
    rho_db = 10; %transmit power [dB]
    angle2 = deg2rad(40); % angle of UE2 to BS
    nbrOfRealizations=1e3; %number of saddlepoint realizations
    Mlist = [2 5 8 10:10:100, 2e2:1e2:5e2 1000 5000]; %the number of antennas considered
    COMBINER = 'M-MMSE'; %what combiner to use [MR, RZF, M-MMSE]
    ESTIMATOR = 'MMSE'; %what combiner to use [LS, MMSE]
    PILOT_CONTAMINATION = 1; %Let the two UEs use the same pilot sequence
    NO_OF_UEs = 2; % 1 or 2
    UNCORRELATED = 0; %0 = spatial correlation
    BS_CSI = 0;
    UE_CSI = 0; %0 = UE uses E[H1'W1], 1 = UE knows H1'W1, 2= UE knows W1 and Hhat1,
    antennaType = 'ULA';
end
% Initializations:
eps_target = 1e-5;
angle1 = deg2rad(30); % angle of UE1 to BS
ASDdeg = 25; %angular delay spread

n_ul = (n-np)/2;
n_dl = n_ul;
rho = 10^(rho_db/10);  %transmit power [linear]

format long
avg_error = inf(length(Mlist), 2);
s_val_ul = nan(length(Mlist));
s_val_dl=nan( length(Mlist));

for i=1:length(Mlist)
     
    M=Mlist(i); %consider current M
    
    if i > 1 % We will evaluate the error probabilities until we achieve the target
        EVALUATE_UL = (avg_error(i-1,1) > eps_target) && ~isinf(avg_error(i-1,1));
        EVALUATE_DL = (avg_error(i-1,2) > eps_target) && ~isinf(avg_error(i-1,2));
    else
        EVALUATE_UL = 1;
        EVALUATE_DL = 1;
    end
    if EVALUATE_DL == 0 && EVALUATE_UL == 0
        continue;
    end
    
    %--------------------------------------------------------
    %Create channel matrices and generate covariance matrices
    %--------------------------------------------------------
    [H1, H2, R1, R2] = generateChannel(M, angle1, angle2, ASDdeg, nbrOfRealizations,UNCORRELATED, antennaType);
    %--------------------------------------------------------------
    % Estimate the channel for UE 1 and UE 2 to the BS serving UE 1
    %--------------------------------------------------------------
    if BS_CSI == 0
        [hhat1, hhat2, C1, C2] =  estimateChannel(rho, np, H1, H2, R1, R2, PILOT_CONTAMINATION,ESTIMATOR);
    else
        hhat1 = H1; C1 = zeros(M,M); hhat2 = H2; C2 = zeros(M,M); %Perfect BS CSI
    end
    %---------------------------------
    %   Create Combining vector
    %---------------------------------
    disp('Combiner')
    tic
    [v1, v2] = createCombiners(rho, hhat1, hhat2, C1, C2, COMBINER, NO_OF_UEs);
    toc
    %--------------------------------------
    % Estimate the uplink error probability
    %--------------------------------------
    tic
    if EVALUATE_UL == 1 
        %--------
        % UPLINK
        %--------
        % Obtain effective noise variance, channel and channel estimation:
        sigma_sq_list = nan(1,nbrOfRealizations);
        g_list = nan(1,nbrOfRealizations);
        g_hat_list= nan(1,nbrOfRealizations);
        for j = 1:nbrOfRealizations
            if NO_OF_UEs == 1
                sigma_sq_list(j) = v1(:,j)'*v1(:,j);
            elseif NO_OF_UEs==2
                sigma_sq_list(j) = rho*abs(v1(:,j)' * H2(:,j))^2 + v1(:,j)'*v1(:,j); 
            end
            g_list(j) = v1(:,j)' * H1(:,j) ;
            g_hat_list(j) = v1(:,j)' * hhat1(:,j);
        end
        % To save time when optimizing s:
        if i == 1
            s_start=1;
        else
            s_start = s_val_ul(i-1); 
        end
        % Computation of the error probability via the saddlepoint approximation:
        [avg_error(i,1), s_val_ul(i)] = sec3_fig2_getErrorProbability(n_ul, rho, b, g_list, g_hat_list, sigma_sq_list, nbrOfRealizations, s_start);
        disp(['UL error probability for (SNR, M) = (' num2str(rho_db) ',' num2str(M) '): ' num2str(avg_error(i,1)) ' and s = ' num2str(s_val_ul(i))]);
    end
    %-----------
    % DOWNLINK
    %-----------
    if EVALUATE_DL == 1
        rate = b / n_dl;
        w1 = v1; %precoding vector for UE1
        w2 = v2; %precoding vector for UE2
        % Obtain effective noise variance, channel and channel estimation:
        if UE_CSI == 0
            ghat_dl = mean(sum(conj(H1).*w1,1)); %UE uses channel hardening
            ghat_dl_list = ghat_dl*ones(1,nbrOfRealizations);
        elseif UE_CSI == 1
            ghat_dl_list = sum(conj(H1).*w1,1); %UE knows the precoded channel
        elseif UE_CSI == 2
            ghat_dl_list = sum(conj(hhat1).*w1,1); %UE knows Hhat1 and assumes it to be correct
        end
        for j = 1:nbrOfRealizations
            if NO_OF_UEs == 1
                sigma_sq_list(j) =  1;
            elseif NO_OF_UEs==2
                sigma_sq_list(j) = rho*abs(H1(:,j)'*w2(:,j))^2 + 1;
            end
            g_list(j) = H1(:,j)'*w1(:,j);
        end
        % To save time when optimizing s:        
        if i == 1
            s_start=1;
        else
            s_start = s_val_dl(i-1) ;
        end
        % Computation of the error probability via the saddlepoint approximation:
        [avg_error(i,2),  s_val_dl(i)] = sec3_fig2_getErrorProbability(n_dl, rho, b, g_list, ghat_dl_list, sigma_sq_list, nbrOfRealizations, s_start);
        disp(['DL error probability for (SNR, M) = (' num2str(rho_db) ',' num2str(M) '): ' num2str(avg_error(i,2)) ' and s = ' num2str(s_val_dl(i))]);
    end
    toc
end

if DEBUG == 1
    
    figure(1); hold on;
    plot(Mlist, avg_error(:,1));
    set(gca, 'YScale', 'log')
    set(gca, 'XScale', 'log')
    ylim([1e-5,1])
    title('UL')
    ylabel('Error probability')
    xlabel('Number of antennas')
    figure (2); hold on;
    plot(Mlist, avg_error(:,2));
    set(gca, 'YScale', 'log')
    set(gca, 'XScale', 'log')
    ylim([1e-5,1])
    title('DL')
    ylabel('Error probability')
    xlabel('Number of antennas')
    
end
avg_error(isinf(avg_error)) = nan;

%--------------------------------------------------------------
%   Save file
%--------------------------------------------------------------
data.np=np;
data.n = n;
data.rate =  rate;
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
data.avg_error_ul = avg_error(:,1);
data.s_val_dl = s_val_dl;
data.avg_error_dl = avg_error(:,2);

filename = [ESTIMATOR '_' COMBINER '_UEs_' num2str(NO_OF_UEs) '_angle2_' num2str(rad2deg(angle2)) '_np_' num2str(np) '_n_' num2str(n) '_ASD_' num2str(ASDdeg) '.mat'];


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
filename = [antennaType '_' filename];
save(filename, 'data', '-v7.3');


end



function [eps_val, s] = asymptoticErrorProbability(angle1, angle2, ASDdeg, antennaSpacing, rho, n, np,rate)
% Function that computes the asymptotic limit as M tends to infinity on the 
% error probability when using MR combining with pilot contamination.
%STEP 1: Find the asymptotic SNR
M=500; % This value turns out to be enough to obtain the asymptotic received SNR
eps_val = ones(size(angle2));
sigma_ul = 1;
n_ul = n-np;

R1 = functionRlocalscattering(M,angle1,ASDdeg,antennaSpacing); %Covariance matrix of UE 1
R2 = functionRlocalscattering(M,angle2,ASDdeg,antennaSpacing); %Covariance matrix of UE 2
Q = (rho*np*R1 + rho*np*R2 + sigma_ul*eye(M)); 
R1Qinv = R1 / Q;
R2Qinv = R2 / Q;
rho_asymp = (real(trace(R1Qinv*R1)/trace(R2Qinv*R1)))^2;
%STEP 2: Find out the smallest epsilon achievable by sweeping s
old_eps_val = inf;
s = 0;
while eps_val <= old_eps_val
    s=s+0.001;
    % Scenario fixed such that the following values hold:
    sigma_sq = 1;
    g = 1;
    ghat = 1;
    % Parameters to compute the saddlepoint approximation:
    betaA_ul = s*rho_asymp*abs(g-ghat)^2 + s*sigma_sq;
    betaB_ul = s*(rho_asymp*abs(g)^2 + sigma_sq) / (1+s*rho_asymp*abs(ghat)^2);
    sigma_v = abs(g)^2 *rho_asymp + sigma_sq;
    gamma = s/(1 + s*rho_asymp*abs(ghat)^2);
    nu_ul = s*gamma*abs(sigma_v - rho_asymp* g'*ghat)^2 / (betaA_ul*betaB_ul);
    preterm_ul = log(1+s*rho_asymp*abs(ghat)^2);
    % Update error probability value and compute the saddlepoint approx:
    old_eps_val = eps_val;
%     [eps_val,~] = saddlepoint_approximation_Alex(n, n_ul, rate, betaA_ul, betaB_ul, nu_ul, preterm_ul);
    eps_val = saddlepoint_approximation(n_ul, rate, betaA_ul, betaB_ul, nu_ul, preterm_ul);
    if s > 100
        break;
    end
end

end

function [H1, H2, R1, R2] = generateChannel(M, angle1, angle2, ASDdeg, nbrOfRealizations, uncorr,antennaType)
% Function that generates the channels and the covariance matrices 
%---------------------------------
% Obtain covariance matrices
%---------------------------------
if uncorr == 0
    if strcmp(antennaType, 'ULA') % Uniform Linear Array
        R1 = functionRlocalscattering(M,angle1,ASDdeg,1/2); %Covariance matrix of UE 1
        R2 = functionRlocalscattering(M,angle2,ASDdeg,1/2); %Covariance matrix of UE 2
        
    elseif strcmp(antennaType, 'UCA') % Uniform Circular Array
        R1 = functionRlocalscatteringApproxUCA(M,angle1,ASDdeg);
        R2 = functionRlocalscatteringApproxUCA(M,angle2,ASDdeg);
    end
    
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
% Function that computes the channel estimations [LS, MMSE]
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
    yp1 = rho*np*H1 + sqrt(np*rho)*Np; %UE 1 and UE2 does not use the same pilot sequence
    yp2 = rho*np*H2 + sqrt(np*rho)*Np; %UE 1 and UE2 does not use the same pilot sequence
    Q1 = (rho*np*R1 + sigma_ul*eye(M));  % matrix for MMSE estimation of H1
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
% Function that computes the combiners and precoders [MR,RZF,M-MMSE]
vecnorm = @(A)  sqrt(sum(A.*conj(A),1)); %norm of each column in matrix

sigma_ul = 1;
M = size(hhat1,1);
nbrOfRealizations = size(hhat1,2);
%---------------------------------
%   Create Combining vector
%---------------------------------
if strcmp(COMBINER, 'MR')
    v1 = hhat1 ./ (vecnorm(hhat1) ); %MR combiner for UE 1
    v2 = hhat2 ./ (vecnorm(hhat2) ); %MR combiner for UE 2
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
