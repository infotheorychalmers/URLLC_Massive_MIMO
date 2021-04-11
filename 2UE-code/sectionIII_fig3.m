function sectionIII_fig3(n, np, rho_db, b, M, nbrOfRealizations, COMBINER, PILOT_CONTAMINATION, UNCORRELATED, NO_OF_UEs, BS_CSI, UE_CSI, nbrOfPositions,antennaType, simNO)
% Function sectionIII_fig2_RCUs_saddlepoint(n, np, rho_db, b, angle2, M, nbrOfRealizations, ESTIMATOR, COMBINER, PILOT_CONTAMINATION, UNCORRELATED, NO_OF_UEs, BS_CSI, UE_CSI): 
% Generates saddlepoint approximations in Fig. 2 for both UL and DL. 
% 
% INPUTS:
% n = blocklength
% np = length pilot sequence
% rho_db = transmit power [dB]
% b = information bits
% M = number of antennas at BS
% nbrOfRealizations = number of realization in the Monte-Carlo simulations
% COMBINER = Combiner to use [MR, M-MMSE, RZF]
% PILOT_CONTAMINATION = [0,1]; 1: the two UEs use the same pilot sequence
% UNCORRELATED = [0,1]; 0: spatially correlated channel
% NO_OF_UEs = number of UEs [1,2]
% BS_CSI = [0,1]; 1: Perfect CSI at BS
% UE_CSI = [0,1,2]; 0: channel hardening; 1: Perfect CSI at UE; 2: UE knows
%                   precoder and channel estimates (same info as in the UL at the BS)
% nbrOfPositions = Number of random positions to generate the CDF
% antennaType = ['ULA','UCA'] uniform linear array or uniform circular array
% simNO: Simulation number (optional). Used externally to generate more batches and get more points of the CDF
DEBUG = 1;

if DEBUG == 1
    %FLAGS
    np = 2; %number of pilots
    n = 300; %total number of channle uses
    b = 8*20; %transmission rate
    rho_db = 10; %transmit power [dBm]
    nbrOfRealizations=1e5; %number of saddlepoint realizations
    nbrOfPositions = 500;
    M = 100; %the number of antennas considered
    COMBINER = 'M-MMSE'; %what combiner to use [MR, M-MMSE]
    PILOT_CONTAMINATION = 1; %Let the two UEs use the same pilot sequence
    NO_OF_UEs = 2; % 1 or 2
    UNCORRELATED = 0; %0 = spatial correlation
    BS_CSI = 0;
    UE_CSI = 0; %0 = UE uses E[H1'W1], 1 = UE knows H1'W1, 2 = UE knows W1 and Hhat1,
    antennaType = 'ULA';
    simNO = 1;
end

% Initializations:
ASDdeg = 25; %angular delay spread
rho = 10^(rho_db/10);  %transmit power [mW]

n_ul = (n-np)/2;
n_dl = n_ul;

L=1; %single cell
K=2; %2 UEs

B = 20e6;%Communication bandwidth
noiseFigure = 7;%Noise figure at the BS (in dB)
noiseVariancedBm = -174 + 10*log10(B) + noiseFigure;%Compute noise power

format long
s_val_ul = nan(1,nbrOfPositions);
s_val_dl=nan(1,nbrOfPositions);
eps_ul = inf(1,nbrOfPositions);
eps_dl = inf(1,nbrOfPositions);

parfor pos_idx = 1:nbrOfPositions
    
    % Initializations to avoid potential issues with the parfor:
    sigma_sq_ul = nan(1,nbrOfRealizations);
    ghat_dl = nan(1,nbrOfRealizations);
    sigma_sq_dl = nan(1,nbrOfRealizations);
    %---------------------------------
    % Obtain covariance matrices
    %---------------------------------
    tic
    [R,channelGaindB] = functionExampleSetup(L,K,M,ASDdeg, antennaType);%Compute channel statistics for one setup 
    channelGainOverNoise = channelGaindB - noiseVariancedBm;%Compute the normalized average channel gain
    
    beta1 = 10^(channelGainOverNoise(1)/10); % Channel gain for UE1 in linear scale
    R1 = beta1 * R(:,:,1); % Correlation matrix UE1 with channel gain
    
    beta2 = 10^(channelGainOverNoise(2)/10); % Channel gain for UE2 in linear scale
    R2 = beta2 * R(:,:,2); % Correlation matrix UE2 with channel gain
    
    if UNCORRELATED == 1
        R1 = beta1*eye(M);
        R2 = beta2*eye(M);
    end
    
    %Small scale fading:
    H1 = sqrt(0.5)*sqrtm(R1)*(randn(M,nbrOfRealizations)+1i*randn(M,nbrOfRealizations));
    H2 = sqrt(0.5)*sqrtm(R2)*(randn(M,nbrOfRealizations)+1i*randn(M,nbrOfRealizations));
    
    % Channel estimation:
    if BS_CSI == 0
        [hhat1, hhat2, C1, C2] =  estimateChannel(rho, np, H1, H2, R1, R2, PILOT_CONTAMINATION);
        
    else
        hhat1 = H1; C1 = zeros(M,M); hhat2 = H2; C2 = zeros(M,M); %Perfect BS CSI
        
    end

    % Compute combiners/precoders
    [v1, v2] = createCombiners(rho, hhat1, hhat2, C1, C2, COMBINER, NO_OF_UEs);
    w1 = v1; %precoding vector for UE1
    w2 = v2; %precoding vector for UE2
    
    % Computation of effective noise, effective channel and channel estimation for the UL:
    g_ul = sum(conj(v1) .* H1,1); % effective channel after combining
    
    ghat_ul = sum(conj(v1) .* hhat1,1); % effective channel estimate after combining
    
    % Compute effective noise:
    if NO_OF_UEs == 1
        sigma_sq_ul = sum(conj(v1) .* v1,1);
    elseif NO_OF_UEs==2
        sigma_sq_ul = rho*abs(sum(conj(v1) .* H2,1)).^2 + sum(conj(v1) .* v1,1);
    end
    
    % Computation of effective noise, effective channel and channel estimation for the DL:
    g_dl = sum(conj(H1) .* w1,1);
    
    if NO_OF_UEs == 1
        sigma_sq_dl =  1;
    elseif NO_OF_UEs==2
        sigma_sq_dl = rho*abs(sum(conj(H1) .* w2,1)).^2 + 1;
    end
    
    if UE_CSI == 0
        ghat_dl = mean(sum(conj(H1).*w1,1)) * ones(1,nbrOfRealizations); %UE uses channel hardening
    elseif UE_CSI == 1
        ghat_dl = sum(conj(H1).*w1,1); %UE knows the precoded channel
    elseif UE_CSI == 2
        ghat_dl = sum(conj(hhat1).*w1,1); %UE knows Hhat1 and assumes it to be correct
    end
    
    time1=toc;
    %-----------------------------------------------
    % Estimate the uplink error probability
    %-----------------------------------------------
    
    tic
    
    rate = b / n_ul; % UL coding rate
    % Computation of the UL average error probability via the saddlepoint
    % approximation:
    f_ul = @(s) getErrorProbabilityUL(s, n_ul, rho, rate, g_ul, ghat_ul, sigma_sq_ul);
    [eps_ul(pos_idx),  s_val_ul(pos_idx)] = searchForCandidateS(f_ul);
    time2 = toc;
    
    %-----------------------------------------------
    % Estimate the downlink error probability
    %-----------------------------------------------
    tic
    rate = b / n_dl; % DL coding rate
    % Computation of the DL average error probability via the saddlepoint
    % approximation:
    f_dl = @(s) getErrorProbabilityDL(s,n_dl, rho, rate, g_dl, ghat_dl, sigma_sq_dl); 
    [eps_dl(pos_idx),  s_val_dl(pos_idx)] = searchForCandidateS(f_dl);
    
    % Print some information about the running times:
    time4 = toc;
    disp(['Combining took: ' num2str(time1) ' seconds'])
    disp(['UL search took: ' num2str(time2) ' seconds'])
    disp(['DL search took: ' num2str(time4) ' seconds'])
    disp(['Position: ' num2str(pos_idx) ' out of ' num2str(nbrOfPositions) ', UL: (eps, s) = (' num2str(eps_ul(pos_idx)) ',' num2str(s_val_ul(pos_idx)) '), DL: (eps, s) = (' num2str(eps_dl(pos_idx)) ',' num2str(s_val_dl(pos_idx)) ')']);
    fprintf('\n')
    
end



if DEBUG == 1
    
    figure(1)
    [cdf, x] = ecdf(eps_ul);
    loglog(x,cdf)
    axis([1e-5 1 0.5 1])
    xlabel('Target error probability')
    ylabel('Network availability')
    figure(2)
    [cdf, x] = ecdf(eps_dl);
    loglog(x,cdf)
    axis([1e-5 1 0.5 1])
    xlabel('Target error probability')
    ylabel('Network availability')
    
    
end

%--------------------------------------------------------------
%   Save file
%--------------------------------------------------------------
data.np=np;
data.n = n;
data.rate =  b/n_ul;
data.snr_db = rho_db;
data.nbrOfRealizations = nbrOfRealizations;
data.nbrOfPositions = nbrOfPositions;
data.M = M;
data.COMBINER = COMBINER;
data.PILOT_CONTAMINATION=PILOT_CONTAMINATION;
data.NO_OF_UEs = NO_OF_UEs;
data.UNCORRELATED = UNCORRELATED;
data.s_val_ul = s_val_ul;
data.avg_error_ul = eps_ul;
data.s_val_dl = s_val_dl;
data.avg_error_dl = eps_dl;

filename = [COMBINER '_M_' num2str(M) '_UEs_' num2str(NO_OF_UEs) '_SNR_' num2str(rho_db) '_np_' num2str(np) '_n_' num2str(n) '.mat'];


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
filename = ['simNr_' num2str(simNO) '_' antennaType '_CDF_' filename];
save(filename, 'data', '-v7.3');


end


function [error_prob, s_val] = searchForCandidateS(f)
% Function to optimize the parameter s of the RCUs
s_list = fliplr(logspace(-8,0,50));

eps_debug = [];

eps_old = inf;
for ii = 1:length(s_list)
    s_candidate = s_list(ii);
    
    eps_cur =f(s_candidate);
    
    eps_debug(ii)=eps_cur;
    
    if eps_cur > eps_old
        eps_cur = eps_old;
        s_candidate = s_list(ii-1);
        break;
    else
        eps_old = eps_cur;
    end
    
end

error_prob = eps_cur;
s_val = s_candidate;
end

function avg_error = getErrorProbabilityUL(s, n, rho,rate, g_list, ghat_list, sigma_sq_list)
% Function getErrorProbabilityUL(s, n, rho, rate, g_list ,ghat_list,sigma_sq_list)
% that computes the quantities needed to compute the CGFs and its
% derivatives to then compute saddlepoiint approximation of the UL average
% error probability given by the RCUs

% Initializations:
nbrOfRealizations = length(g_list);
eps_ul = nan(1, nbrOfRealizations);
for j = 1:nbrOfRealizations
    % Get channel, channel estimate, and effective noise:
    g = g_list(j);
    ghat = ghat_list(j);
    sigma_sq = sigma_sq_list(j);
    % Parameters related to the CGF of the info density:
    betaA_ul = s*rho*abs(g-ghat)^2 + s*sigma_sq;
    betaB_ul = s*(rho*abs(g)^2 + sigma_sq) / (1+s*rho*abs(ghat)^2);
    sigma_v = abs(g)^2 *rho + sigma_sq;
    gamma = s/(1 + s*rho*abs(ghat)^2);
    nu_ul = s*gamma*abs(sigma_v - rho* g'*ghat)^2 / (betaA_ul*betaB_ul);
    preterm_ul = log(1+s*rho * abs(ghat)^2);
    % Compute the saddlepoint approximation
    eps_ul(j) = saddlepoint_approximation(n, rate, betaA_ul, betaB_ul, nu_ul, preterm_ul);   
end
avg_error=mean(eps_ul); % Average over all the random realizations
end

function avg_error = getErrorProbabilityDL(s, n, rho,rate, g_list, ghat_list, sigma_sq_list)
% Function getErrorProbabilityDL(s, n, rho, rate, g_list ,ghat_list,sigma_sq_list)
% that computes the quantities needed to compute the CGFs and its
% derivatives to then compute saddlepoiint approximation of the DL average
% error probability given by the RCUs

% Initializations:
nbrOfRealizations = length(g_list);
%check channel hardening for debugging purposes:
% g = sum(conj(H1).*w1,1);
% ghat = mean(sum(conj(H1).*w1,1));
% plot(abs(g),'b'); hold on;plot(1:length(g),ones(size(g))*abs(ghat),'r--')

eps_dl = nan(1, nbrOfRealizations);
for j = 1:nbrOfRealizations
    % Get channel, channel estimate, and effective noise:
    g = g_list(j);
    ghat = ghat_list(j);
    sigma_sq = sigma_sq_list(j);
    % Parameters related to the CGF of the info density:      
    betaA_dl = s*rho*abs(g-ghat)^2 + s*sigma_sq;
    betaB_dl = s*(rho*abs(g)^2 + sigma_sq) / (1+s*rho*abs(ghat)^2);
    sigma_v = abs(g)^2 *rho + sigma_sq;
    gamma = s/(1 + s*rho*abs(ghat)^2);
    nu_dl = s*gamma*abs(sigma_v - rho* g'*ghat)^2 / (betaA_dl*betaB_dl);
    preterm_dl = log(1+s*rho * abs(ghat)^2);
    % Compute the saddlepoint approximation
    eps_dl(j) = saddlepoint_approximation(n, rate, betaA_dl, betaB_dl, nu_dl, preterm_dl);
    
end
avg_error=  mean(eps_dl); % Average over all the random realizations
end

function [hhat1, hhat2, C1, C2] =  estimateChannel(rho, np, H1, H2, R1, R2, PILOT_CONTAMINATION)
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
    yp1 = rho*np*H1 + sqrt(np*rho)*Np; %UE 1 and UE2 does not use the same pilot sequence
    yp2 = rho*np*H2 + sqrt(np*rho)*Np; %UE 1 and UE2 does not use the same pilot sequence
    Q1 = (rho*np*R1 + sigma_ul*eye(M));  % matrix for MMSE estimation of H1
    Q2 = (rho*np*R2 + sigma_ul*eye(M)); % matrix for MMSE estimation of H2
end
R1Qinv = R1 / Q1;
R2Qinv = R2 / Q2;
C1 = R1 - np*rho*R1Qinv*R1; %covariance of Hhat1
C2 = R2 - np*rho*R2Qinv*R2; %covariance of Hhat2
hhat1 = R1Qinv*yp1; %estimate H1
hhat2 = R2Qinv*yp2; %estimate H2
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
elseif strcmp(COMBINER, 'M-MMSE')   
    v1 = nan(size(hhat1));
    v2 = nan(size(hhat2));
    if NO_OF_UEs == 1
        Z = C1 + sigma_ul*eye(M); %matrix used for the inverse in M-MMSE estimation       
        for k = 1:nbrOfRealizations %estimate each of the realizations
            v1(:,k) = rho* ( (rho*(hhat1(:,k)*hhat1(:,k)') + Z) \ hhat1(:,k)); % M-MMSE combiner for UE 1
            v1(:,k) = v1(:,k)/norm(v1(:,k));
        end        
    elseif NO_OF_UEs == 2
        Z = C1 + C2 + sigma_ul*eye(M); %matrix used for the inverse in M-MMSE estimation        
        for k = 1:nbrOfRealizations %estimate each of the realizations
            PHI = (rho*(hhat1(:,k)*hhat1(:,k)') + rho*(hhat2(:,k)*hhat2(:,k)') + Z);            
            v1(:,k) = rho* ( PHI \ hhat1(:,k)); % M-MMSE combiner for UE 1
            v1(:,k) = v1(:,k)/norm(v1(:,k));            
            v2(:,k) = rho * ( PHI \ hhat2(:,k)); % M-MMSE combiner for UE 2
            v2(:,k) = v2(:,k)/norm(v2(:,k));
        end        
    end
end
end
