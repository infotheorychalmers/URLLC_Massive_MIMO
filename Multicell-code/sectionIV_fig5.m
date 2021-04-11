function data = sectionIV_fig5(L, K, M, infobits, n, np, rho_db, ESTIMATOR, COMBINER, nbrOfRealizations,nbrOfPositions, antennaType, simNO, UE_CSI)
% Function sectionIV_fig5(n, np, rho_db, b, angle2, M, nbrOfRealizations, ESTIMATOR, COMBINER, PILOT_CONTAMINATION, UNCORRELATED, NO_OF_UEs, BS_CSI, UE_CSI): 
% Generates saddlepoint approximations in Fig. 2 for both UL and DL. 
% 
% INPUTS:
% L = number of cells
% K = number of users per cell
% M = number of antennas at BS
% infobits = information bits
% n = blocklength
% np = length pilot sequence
% rho_db = transmit power [dB]
% ESTIMATOR = Channel estimator to use [LS, MMSE]
% COMBINER = Combiner to use [MR, M-MMSE, RZF]
% nbrOfRealizations = number of realization in the Monte-Carlo simulations
% nbrOfPositions = Number of random positions to generate the CDF
% antennaType = ['ULA','UCA'] uniform linear array or uniform circular array
% simNO: Simulation number (optional). Used externally to generate more batches and get more points of the CDF
% UE_CSI = [0,1,2]; 0: channel hardening; otherwise: UE knows precoder and
%                   channel estimates (same info as in the UL at the BS)


vecnorm = @(A)  sqrt(sum(A.*conj(A),1)); %norm of each column in matrix

DEBUG = 1;
rng('shuffle')
if nargin < 14
    UE_CSI=0;
end

if DEBUG == 1
    %FLAGS
    L=4;
    K=10;
    np = 40; %number of pilots
    n = 300; %total number of channle uses
    infobits = 8*20; %transmission rate
    rho_db = 10; %transmit power [dBm]
    nbrOfRealizations=1e2; %number of saddlepoint realizations
    nbrOfPositions = 1;
    M = 100; %the number of antennas considered
    COMBINER = 'M-MMSE'; %what combiner to use [MR, M-MMSE, RZF]
    ESTIMATOR = 'MMSE'; %what combiner to use [LS, MMSE]
    antennaType = 'ULA';
    simNO = 1;
end
% Preliminary initializations: 
eps_target = 1e-6; % We don't continue optimizations if est. epsilon < 1e-6
n_ul = (n - np) / 2; %length of UL phase
n_dl = n - np - n_ul; %length of DL phase


%% Propagation parameters
B = 20e6;%Communication bandwidth
noiseFigure = 7;%Noise figure at the BS (in dB)
noiseVariancedBm = -174 + 10*log10(B) + noiseFigure;%Compute noise power
ASDdeg = 25; %angular delay spread

p = 10^(rho_db/10); % transmission power in linear scale
%% Run Simulation

% Initializations:
eps_ul = inf(nbrOfPositions,1);
eps_dl = inf(nbrOfPositions,1);
s_val_ul = inf(nbrOfPositions,1);
s_val_dl = inf(nbrOfPositions,1);

parfor pos_idx = 1:nbrOfPositions % Go through all setups
    
    %Create batches: (to avoid RAM issues)
    tic
    [R,channelGaindB, posXY, posBS] = functionExampleSetup(L,K,M,ASDdeg,antennaType);%Compute channel statistics for one setup
    channelGainOverNoise = channelGaindB - noiseVariancedBm;%Compute the normalized average channel gain
    
    
    k = randi(K,1); %generate a random UE
    j = randi(L,1); %generate a random cell
    ue_indx = (j-1)*K + k; %which one of all L*K UEs are we considering
    
    % Further initializations:
    batchsize = min(1e2, nbrOfRealizations);
    batches = nbrOfRealizations / batchsize;
    g_ul = nan(nbrOfRealizations,1);
    ghat_ul = g_ul;
    sigma_sq_ul = g_ul;
    g_dl = g_ul;
    ghat_dl = g_ul;
    sigma_sq_dl = g_ul;
    %
    for batch_indx =1:batches %Need to do this for memory issues
        %Generate channel realizations with estimates and estimation error correlation matrices:
        [Hhat,C,Rscaled, H] = functionChannelEstimates(R,channelGainOverNoise,batchsize,M,K,L,p,np, ESTIMATOR);
        %Create combining/precoding vectors:
        [~, W] =  functionCombinerPrecoder(Hhat,C,Rscaled,batchsize,M,K,L,p, COMBINER);
        
        v_jk = W(:,:,k,j);  %combiner vector for UE k in cell j (normalized to make optimal s smaller)
        hhat_jkj = Hhat(:,:,k,j,j); %estimated channel between  BS in cell j and UE k in cell j (hhat_{j,k}^j)
        h_jkj = H(:,:,k,j,j); %channel between UE k in cell j to BS j, size: M x batchsize
        Hj = squeeze(H(:,:,:,:,j)); % all channels between the UEs and BS j, size: M x batchsize x K x L

        for batch_sample = 1:batchsize
            %UL
            Hj_samp=squeeze(Hj(:,batch_sample,:,:)); %pick out current sample, size: M x K x L
            Hj_samp = reshape(Hj_samp, size(Hj_samp,1), size(Hj_samp,2)*size(Hj_samp,3)); %H^j = [h_11^j, ..., h_LK^j], size: M x KL

            combinedChannels = v_jk(:,batch_sample)' * Hj_samp;

            g_ul((batch_indx-1)*batchsize + batch_sample) = combinedChannels(ue_indx); % effective channel after combining
            ghat_ul((batch_indx-1)*batchsize +batch_sample) = v_jk(:,batch_sample)' * hhat_jkj(:,batch_sample); % effective channel estimate after combining
            combinedChannels(ue_indx) = []; %remove the UEs channel to only have noise terms left
            sigma_sq_ul((batch_indx-1)*batchsize + batch_sample) = vecnorm(v_jk(:,batch_sample))^2 +  p*sum(abs(combinedChannels).^2); % effective noise
            
            %DL
            combinedChannels = nan(1, K*L);
            for l = 1:L %go through each cell
                combinedChannels((l-1)*K+1:l*K) = H(:,batch_sample,k,j,l)' * squeeze(W(:,batch_sample,:,l)); % (h_kj^l)^H * [w_1^l, ... , w_K^l]
            end
            
            g_dl((batch_indx-1)*batchsize + batch_sample) = combinedChannels(ue_indx); % effective channel after combining
            combinedChannels(ue_indx) = []; %remove the UEs channel to only have noise terms left
            sigma_sq_dl((batch_indx-1)*batchsize + batch_sample) = 1 +  p*sum(abs(combinedChannels).^2); % effective noise
        end
        ghat_dl((batch_indx-1)*batchsize + 1 : batch_indx*batchsize) = sum(conj(h_jkj) .* W(:,:,k,j), 1); % this will be used to estimate the channel hardening for DL

    end
    
    if UE_CSI == 0
        ghat_dl = mean(ghat_dl) * ones(size(g_dl)); %estimate channel hardening
    end
    
    time = toc;
    disp(['Combiner + estimation takes: ' num2str(time) ' seconds']);

    %-----------------------------------------------
    % Estimate the uplink error probability
    %-----------------------------------------------
    tic   
    rate = infobits / n_ul; % UL coding rate
    f_ul = @(s) getErrorProbabilityUL(s, n_ul, p, rate, g_ul, ghat_ul, sigma_sq_ul); % function of s in the saddlepoint approximation of the error 
    [eps_ul(pos_idx),  s_val_ul(pos_idx)] = searchForCandidateS(f_ul,eps_target); % Optimizaton over s
    time2 = toc;
    disp(['UL search took: ' num2str(time2) ' seconds'])
    %-----------------------------------------------
    % Estimate the downlink error probability
    %-----------------------------------------------
    rate = infobits / n_dl; %DL coding rate
    tic
    f_dl = @(s) getErrorProbabilityDL(s,n_dl, p, rate, g_dl, ghat_dl, sigma_sq_dl); 
    [eps_dl(pos_idx),  s_val_dl(pos_idx)] = searchForCandidateS(f_dl,eps_target);
    time4=toc;
    disp(['DL search took: ' num2str(time4) ' seconds'])
    disp(['UL: epsilon = ' num2str(eps_ul(pos_idx)) ' and s = ' num2str(s_val_ul(pos_idx))]);
    disp(['DL: epsilon = ' num2str(eps_dl(pos_idx)) ' and s = ' num2str(s_val_dl(pos_idx))]);
    fprintf('\n')

end


if DEBUG == 1
    
    subplot(1,2,1)
    [cdf, x]=ecdf(eps_ul);
    loglog(x,cdf); hold on;
    axis([1e-5 1 0.5 1])
    xlabel('Target error probability')
    ylabel('Network availability')
    
    subplot(1,2,2)
    [cdf, x]=ecdf(eps_dl);
    loglog(x,cdf); hold on;
    axis([1e-5 1 0.5 1])
    xlabel('Target error probability')
    ylabel('Network availability')
    
    suptitle(['M=' num2str(M) ', K=', num2str(K), ', tx power=' num2str(rho_db)])
end

%--------------------------------------------------------------
%   Save file
%--------------------------------------------------------------
data.ASDdeg = ASDdeg;
data.L=L;
data.K=K;
data.M=M;
data.np=np;
data.n = n;
data.rate =  infobits/n_ul;
data.snr_db = rho_db;
data.nbrOfRealizations = nbrOfRealizations;
data.nbrOfPositions = nbrOfPositions;
data.Mlist = M;
data.COMBINER = COMBINER;
data.ESTIMATOR = ESTIMATOR;
data.s_val_ul = s_val_ul;
data.avg_error_ul = eps_ul;
data.s_val_dl = s_val_dl;
data.avg_error_dl = eps_dl;

filename = [ESTIMATOR '_' COMBINER '_M_' num2str(M) '_SNR_' num2str(rho_db) '_np_' num2str(np) '_n_' num2str(n) '_ASD_' num2str(ASDdeg) '_K_' num2str(K) '_L_' num2str(L) '.mat'];
filename = ['simNr_' num2str(simNO) '_Multicell_' antennaType '_CDF_' filename];
save(filename, 'data', '-v7.3');

end


function [error_prob, s_val] = searchForCandidateS(f,eps_target)
% Function to optimize the parameter s of the RCUs
s_list = fliplr(logspace(-8,0,50));

% eps_debug = []; % for debugging purposes

eps_old = inf;
for ii = 1:length(s_list)
    s_candidate = s_list(ii);
    
    eps_cur =f(s_candidate);
    
    if eps_cur < eps_target      
        break;
    end
    
    %eps_debug(ii)=eps_cur; % for debugging purposes
    
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

function [avg_error] = getErrorProbabilityUL(s, n, rho, rate, g_list, ghat_list, sigma_sq_list)
% Function getErrorProbabilityUL(s, n, rho, rate, g_list ,ghat_list,sigma_sq_list)
% that computes the quantities needed to compute the CGFs and its
% derivatives to then compute saddlepoiint approximation of the UL average
% error probability given by the RCUs

% Initializations:
nbrOfRealizations = length(g_list);
eps_ul = nan(1, nbrOfRealizations);
for i = 1:nbrOfRealizations
    % Get channel, channel estimate, and effective noise:
    g = g_list(i);
    ghat = ghat_list(i);
    sigma_sq= sigma_sq_list(i);
    % Parameters related to the CGF of the info density:
    betaA_ul = s*rho*abs(g-ghat)^2 + s*sigma_sq;
    betaB_ul = s*(rho*abs(g)^2 + sigma_sq) / (1+s*rho*abs(ghat)^2);
    sigma_v = abs(g)^2 *rho + sigma_sq;
    gamma = s/(1 + s*rho*abs(ghat)^2);
    nu_ul = s*gamma*abs(sigma_v - rho* g'*ghat)^2 / (betaA_ul*betaB_ul);
    preterm_ul = log(1+s*rho * abs(ghat)^2);
    % Compute the saddlepoint approximation
    eps_ul(i) = saddlepoint_approximation(n, rate, betaA_ul, betaB_ul, nu_ul, preterm_ul);   
end
avg_error=mean(eps_ul); % Average over all the random realizations

end

function [avg_error] = getErrorProbabilityDL(s, n, rho, rate, g_list, ghat_list, sigma_sq_list)
% Function getErrorProbabilityDL(s, n, rho, rate, g_list ,ghat_list,sigma_sq_list)
% that computes the quantities needed to compute the CGFs and its
% derivatives to then compute saddlepoiint approximation of the DL average
% error probability given by the RCUs

% Initializations:
nbrOfRealizations = length(g_list);
%check channel hardening
%  g_test = sum(conj(H(:,:,k,j,j)) .* squeeze(W(:,:,:,j)),1);
%  plot(abs(g_test),'b'); hold on;plot(1:length(g_test),ones(size(g_test))*abs(g_hat_dl),'r--')

eps_dl = nan(1, nbrOfRealizations);
for i = 1:nbrOfRealizations
    % Get channel, channel estimate, and effective noise:
    g = g_list(i);
    ghat = ghat_list(i); 
    sigma_sq= sigma_sq_list(i);
    % Parameters related to the CGF of the info density:
    betaA_dl = s*rho*abs(g-ghat)^2 + s*sigma_sq;
    betaB_dl = s*(rho*abs(g)^2 + sigma_sq) / (1+s*rho*abs(ghat)^2);
    sigma_v = abs(g)^2 *rho + sigma_sq;
    gamma = s/(1 + s*rho*abs(ghat)^2);
    nu_dl = s*gamma*abs(sigma_v - rho* g'*ghat)^2 / (betaA_dl*betaB_dl);
    preterm_dl = log(1+s*rho * abs(ghat)^2);
    % Compute the saddlepoint approximation
    eps_dl(i) = saddlepoint_approximation(n, rate, betaA_dl, betaB_dl, nu_dl, preterm_dl);    
end
avg_error=  mean(eps_dl); % Average over all the random realizations
end
