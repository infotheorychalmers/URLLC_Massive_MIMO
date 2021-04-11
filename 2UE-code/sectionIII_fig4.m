function sectionIII_fig4(n, np, rho_db, rate, angle2, Mlist, nbrOfRealizations, COMBINER, PILOT_CONTAMINATION, PERFECT_CSI, NO_OF_UEs, OUTAGE)
% Function sectionIII_fig4(n, np, rho_db, rate, angle2, Mlist, nbrOfRealizations, COMBINER, PILOT_CONTAMINATION, PERFECT_CSI, NO_OF_UEs, OUTAGE): 
% Generates saddlepoint approximations in Fig. 2 for both UL and DL. 
% 
% INPUTS:
% n = blocklength
% np = length pilot sequence
% rho_db = transmit power [dB]
% rate = rate in bits
% Mlist = list with number of antennas at BS
% nbrOfRealizations = number of realization in the Monte-Carlo simulations
% COMBINER = Combiner to use [MR, M-MMSE, RZF]
% PILOT_CONTAMINATION = [0,1]; 1: the two UEs use the same pilot sequence
% PERFECT_CSI = [0,1]; 1: perfect CSI (only valid for MR)
% NO_OF_UEs = number of UEs [1,2]
% OUTAGE = [0,1]; 1: Evaluate outage error probability instead of RCUs

DEBUG = 1;

if DEBUG == 1
    %FLAGS
    np = 2; %number of pilots
    n = 300; %total number of channle uses
    n_ul =((n-np)/2);
    rate = 20*8/n_ul; %transmission rate
    rho_db = 10; %transmit power [dB]
    angle2 = deg2rad(40); % angle of UE2 to BS
    nbrOfRealizations=1e4; %number of saddlepoint realizations
    Mlist = [1 5 7 8]; %the number of antennas considered
    COMBINER = 'M-MMSE'; %what combiner to use [MR, M-MMSE]
    PILOT_CONTAMINATION = 0; %Let the two UEs use the same pilot sequence
    PERFECT_CSI = 0; %perfect CSI (only for MR)
    NO_OF_UEs = 2; % 1 or 2
    OUTAGE = 0; %evaluate outage error probability instead of RCUs
end
if COMBINER == 0
    COMBINER = 'MR';
elseif COMBINER == 1
    COMBINER = 'M-MMSE';
end

vecnorm = @(A)  sqrt(sum(A.*conj(A),1)); %norm of each column in matrix

rho = 10^(rho_db/10);  %transmit power [linear]
angle1 = deg2rad(30); % angle of UE1 to BS

ASDdeg = 25; %angular delay spread
antennaSpacing = 1/2; %antenna separation
sigma_ul = 1; %additive noise variance

if OUTAGE == 1
    nbrOfRealizations=1e7;
    PILOT_CONTAMINATION = 0;
    PERFECT_CSI = 1;
end

if NO_OF_UEs == 1
    PILOT_CONTAMINATION = 0; %only one UE so no pilot contamination
end

format long
avg_error = inf(length(Mlist), 2);

FLAG_UL = 1;
FLAG_DL = 1;

s_val_ul = nan(size(Mlist));
s_val_dl=nan(size(Mlist));

% For each value of M:
for i=1:length(Mlist)
    M=Mlist(i); %consider current M
    
    %---------------------------------
    % Obtain covariance matrices
    %---------------------------------
    R1 = functionRlocalscattering(M,angle1,ASDdeg,antennaSpacing); %Covariance matrix of UE 1
    R2 = functionRlocalscattering(M,angle2,ASDdeg,antennaSpacing); %Covariance matrix of UE 2
    
    %---------------------------------
    %Create channel matrices
    %---------------------------------
    H1 = sqrt(0.5)*sqrtm(R1)*(randn(M,nbrOfRealizations)+1i*randn(M,nbrOfRealizations));
    H2 = sqrt(0.5)*sqrtm(R2)*(randn(M,nbrOfRealizations)+1i*randn(M,nbrOfRealizations));
    
    %--------------------------------------------------------------
    % Estimate the channel for UE 1 and UE 2 to the BS serving UE 1
    %--------------------------------------------------------------
    if PERFECT_CSI == 1
        PILOT_CONTAMINATION = 0;
        hhat1 = H1; hhat2 = H2;
        np = 0; n_ul=(n-np)/2;
        C1=0; C2=0;
    else
        n_ul = (n-np)/2; %number of channel uses left for data
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
    %--------------------------------------
    % Estimate the uplink error probability
    %--------------------------------------
    if FLAG_UL == 1
        FLAG_UL = (i==1) || (avg_error(i-1,1) > 1e-9);
    end
    if FLAG_UL == 1 %only interested in epsilon > 1e-9
        if OUTAGE == 1
            avg_error(i) = getOutageErrorProbability(H1, H2, v1, rho, rate, NO_OF_UEs);
            s_val_ul = nan;
        else

            f = @(s) getErrorProbabilityUL(s, n_ul, rho, rate, hhat1, H1, H2, v1, NO_OF_UEs);
            if i > 1 
                END_INT = s_val_ul(i-1);
                START_INT = END_INT/100;
                TOL = START_INT;
            else
                END_INT = 1;
                START_INT = 1e-3;
                TOL = 1e-4;
            end

            [avg_error(i,1), s_val_ul(i)] = golden_search(f, START_INT, END_INT, TOL);

            disp(['UL error probability for M = ' num2str(M) ': ' num2str(avg_error(i,1)) ' and s = ' num2str(s_val_ul(i))]);
        end
    end
    
    %----------------------------------------
    % Estimate the downlink error probability
    %----------------------------------------
    if FLAG_DL == 1
        FLAG_DL = (i==1) || (avg_error(i-1,2) > 1e-9);
    end
    
    if FLAG_DL == 1 %only interested in epsilon > 1e-9

        f = @(s) getErrorProbabilityDL(s, n_ul, rho,rate, hhat1, H1, v1, v2, NO_OF_UEs);
        if i > 1
            END_INT = s_val_dl(i-1);
            START_INT = END_INT/100;
            TOL = START_INT;
        else
            END_INT = 1;
            START_INT = 1e-3;
            TOL = 1e-4;
        end
        
        [avg_error(i,2), s_val_dl(i)] = golden_search(f, START_INT, END_INT, TOL);
        
        disp(['DL error probability for M = ' num2str(M) ': ' num2str(avg_error(i,2)) ' and s = ' num2str(s_val_dl(i))]);
        
    end
    
    
end
if DEBUG == 1
    keyboard;
    figure(1)
    loglog(Mlist, avg_error,'-*'); hold on;
    axis([1 500 1e-5 1])
    xlabel('number of antennas')
    ylabel('error probability')
end
avg_error(isinf(avg_error)) = nan;
%--------------------------------
%   MR Pilot Contamination limit:
%--------------------------------
if PILOT_CONTAMINATION == 1 && strcmp(COMBINER, 'MR') == 1
    error_prob_asymptotic = asymptoticErrorProbability(angle1, angle2, ASDdeg, antennaSpacing, rho,n, np, rate);
    
    if DEBUG == 1
        if length(angle2) > 1
            figure(2)
            semilogy(rad2deg(angle2), error_prob_asymptotic);
            xlabel('Angle of UE 2')
        end
        figure(1)
        keyboard
        loglog(Mlist, error_prob_asymptotic*ones(size(Mlist)), '*--')
        %ylim([error_prob_asymptotic/2,1])
        title(['SNR=' num2str(rho_db) ', UE2 angle = ' num2str(rad2deg(angle2))])
    end
    data.error_prob_asymptotic=error_prob_asymptotic;
end
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
data.PILOT_CONTAMINATION=PILOT_CONTAMINATION;
data.PERFECT_CSI=PERFECT_CSI;
data.NO_OF_UEs = NO_OF_UEs;
data.OUTAGE = OUTAGE;
data.s_val_ul = s_val_ul;
data.avg_error_ul = avg_error(:,1);
data.s_val_dl = s_val_dl;
data.avg_error_dl = avg_error(:,2);

filename = [COMBINER '_UEs_' num2str(NO_OF_UEs) '_SNR_' num2str(rho_db) '_angle2_' num2str(rad2deg(angle2)) '_np_' num2str(np) '_n_' num2str(n) '.mat'];
if OUTAGE == 1 
    filename = ['Outage_' filename];
elseif PERFECT_CSI == 1
    filename = ['CSI_' filename];
end

if PERFECT_CSI == 0
    filename = ['imperfect_CSI_' filename];
end

if PILOT_CONTAMINATION == 1
    filename = ['pilot_contamination_' filename];
end

save(['sectionIII_fig4_' filename], 'data', '-v7.3');


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

function avg_error = getErrorProbabilityUL(sul, n_ul, rho,rate, hhat1, H1, H2, v1, NO_OF_UEs)
% Function that computes the average error probability dor the UL

% Initializations:
nbrOfRealizations = size(hhat1,2);
eps_ul = nan(1, nbrOfRealizations);

parfor j = 1:nbrOfRealizations
    sigma_sq = 1;
    if NO_OF_UEs == 1
        sigma_sq = 1;%v1 must be normalized to 1 for this to hold.
    elseif NO_OF_UEs==2
        sigma_sq = rho*abs(v1(:,j)' * H2(:,j))^2 + 1; %v1 must be normalized to 1 for this to hold.
    end
    % Obtain effective channel and channel estimation after combining   
    g = v1(:,j)' * H1(:,j);
    ghat = v1(:,j)' * hhat1(:,j);
    % Parameters for the saddlepoint approximation:
    betaA_ul = sul*rho*abs(g-ghat)^2 + sul*sigma_sq;
    betaB_ul = sul*(rho*abs(g)^2 + sigma_sq) / (1+sul*rho*abs(ghat)^2);
    sigma_v = abs(g)^2 *rho + sigma_sq;
    gamma = sul/(1 + sul*rho*abs(ghat)^2);
    nu_ul = sul*gamma*abs(sigma_v - rho* g'*ghat)^2 / (betaA_ul*betaB_ul);
    preterm_ul = log(1+sul*rho * abs(ghat)^2);
    % Compute the error probability via the saddlepoint approximation:
%     [eps_ul(j), reg_ul(j) ] = saddlepoint_approximation_Alex(n, n_ul, rate, betaA_ul, betaB_ul, nu_ul, preterm_ul);
    eps_ul(j) = saddlepoint_approximation(n_ul, rate, betaA_ul, betaB_ul, nu_ul, preterm_ul);
    
end
avg_error=mean(eps_ul);

end

function avg_error = getErrorProbabilityDL(sul, n_dl, rho,rate, hhat1, H1, v1, v2, NO_OF_UEs)
% Function that computes the average error probability dor the UL

% Initializations:
nbrOfRealizations = size(hhat1,2);
eps_dl = nan(1, nbrOfRealizations);

% Precoders equal to the combiners due to UL-DL duality
w1 = v1;
w2 = v2;
ghat_dl = mean(sum(conj(H1).*w1,1)); % channel estimate using channel hardening

parfor j = 1:nbrOfRealizations
    sigma_sq = 1;
    if NO_OF_UEs == 1
        sigma_sq = 1;
    elseif NO_OF_UEs==2
        sigma_sq = rho*abs(H1(:,j)'*w2(:,j))^2 + 1;
    end
    % Obtain effective channel and channel estimation after combining 
    g = H1(:,j)'*w1(:,j);
    ghat = ghat_dl;
    % Parameters for the saddlepoint approximation:
    betaA_dl = sul*rho*abs(g-ghat)^2 + sul*sigma_sq;
    betaB_dl = sul*(rho*abs(g)^2 + sigma_sq) / (1+sul*rho*abs(ghat)^2);
    sigma_v = abs(g)^2 *rho + sigma_sq;
    gamma = sul/(1 + sul*rho*abs(ghat)^2);
    nu_dl = sul*gamma*abs(sigma_v - rho* g'*ghat)^2 / (betaA_dl*betaB_dl);
    preterm_dl = log(1+sul*rho * abs(ghat)^2);
    % Compute the error probability via the saddlepoint approximation:
%     [eps_dl(j), reg_dl(j) ] = saddlepoint_approximation_Alex(n, n_ul, rate, betaA_dl, betaB_dl, nu_dl, preterm_dl);
    eps_dl(j) = saddlepoint_approximation(n_dl, rate, betaA_dl, betaB_dl, nu_dl, preterm_dl);
end
avg_error=  mean(eps_dl);

end

