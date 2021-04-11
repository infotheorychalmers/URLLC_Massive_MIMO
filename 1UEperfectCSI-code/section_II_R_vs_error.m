function section_II_R_vs_error(n,b_vec,snr_db,M,nbrOfRealizations)
% Function section_II_R_vs_error(n,b_vec,snr_db,M,nbrOfRealizations): Generates 
% a rate vs error probability plot for a given number of antennas.
% 
% INPUTS:
% n = blocklength
% b_vec = list with values of information bits
% snr_db = SNR [dB]
% Mlist = number of antennas
% nbrOfRealizations = number of realization in the Monte-Carlo simulations
%
DEBUG = 1;
if DEBUG == 1
    snr_db = -24;
    M = 200;
    nbrOfRealizations = 1e4;
    n = 100;
    b_vec = [20:10:140];
end
R_vec = b_vec/n;
rho = 10.^(snr_db/10);%rho = avergage gain * transmit power / sigma^2
% Definition of functions:
vecnorm = @(A)  sqrt(sum(A.*conj(A),1)); %norm of each column in matrix
% Initializations:
eps_out = nan(1,length(R_vec));
eps_rcus_sp = nan(1,length(R_vec));
eps_rcus = nan(1,length(R_vec));
eps_na = nan(1,length(R_vec));
I1 = nan;
% For golden search:
START_INT = 1e-2;
END_INT = 1;
TOL = 1e-2;
%-------------------------------------
% GENERATE CHANNEL
H1 = generateChannel(M, nbrOfRealizations);
hhat1 = H1;
v1 = hhat1 ./ (vecnorm(hhat1)); %MR combiner for UE 1hhat1 = H1;
sigma_sq = nan(1,nbrOfRealizations);
g = nan(1,nbrOfRealizations);
ghat = nan(1,nbrOfRealizations);

for j=1:nbrOfRealizations
    sigma_sq(j) = v1(:,j)'*v1(:,j);
    g(j) = v1(:,j)' * H1(:,j) ;
    ghat(j) = v1(:,j)' * hhat1(:,j);
end
%-------------------------------------
%COMPUTE Is for chosen s: For Perfect CSI s=1 optimizes GMI.
betaA_ul = rho*abs(g-ghat).^2 + sigma_sq;
betaB_ul = (rho*abs(g).^2 + sigma_sq) ./ (1+rho*abs(ghat).^2);
preterm_ul = log(1+rho * abs(ghat).^2);
I1 = mean(preterm_ul + betaA_ul - betaB_ul)/log(2);
for i = 1:length(R_vec)
    R = R_vec(i);
    b = b_vec(i);
    disp(['R = ' num2str(R)])
    %-------------------------------------
    % EVALUATE RCUs for perfect CSI
    f = @(s) getErrorProbability(s, n, rho, R, g,ghat,sigma_sq);
    [eps_rcus_sp(i), s_val] = golden_search(f, START_INT, END_INT, TOL);
    if eps_rcus_sp(i) > 1e-3 % We only compute real RCUs for large epsilons
        [nz,nv,ng]=generateRVs(rho, n, M, nbrOfRealizations);
        i_s = @(s) -s*nz + (s./(1+s*rho*ng)).*nv + n*log(1+s*rho*ng);
        error_val = @(s) mean(exp(-max(0, i_s(s) - log(2^b-1))));
        [eps_rcus(i), s_val(i)]=golden_search(error_val,  1e-2, 2, 1e-4);
    end
        
    %If we want to generate points of the exact RCUs independently of
    %the results of the saddlepoint:
    
    % [nz,nv,ng]=generateRVs(rho, n, M, nbrOfRealizations);
    % i_s = @(s) -s*nz + (s./(1+s*rho*ng)).*nv + n*log(1+s*rho*ng);
    % error_val = @(s) mean(exp(-max(0, i_s(s) - log(2^b-1))));
    % [eps_rcus(i), s_val(i)]=golden_search(error_val,  1e-2, 2, 1e-4);
    
    %-------------------------------------
    %EVALUATE OUTAGE
    eps_out(i) = gamcdf((2^R - 1)/rho, M, 1); %compute outage probability
    %-------------------------------------
    %EVALUATE NORMAL APPROXIMATION
    fNA = @(s) getErrorProbabilityNA(s, n, rho, b, g,ghat,sigma_sq); % Arguments for numerical option: (s, n, rho, b, H1)
    [eps_na(i), s_val_na(i)]=golden_search(fNA,  START_INT, END_INT, TOL);    
end
%-------------------------------------
% FIGURES
%-------------------------------------
if DEBUG == 1
figure
semilogy(R_vec, eps_na); hold on;
semilogy(R_vec, eps_rcus,'black*');
semilogy(R_vec, eps_rcus_sp,'black');
semilogy(R_vec, eps_out,'r');
semilogy(ones(size(R_vec))*I1, linspace(1e-7,1,length(R_vec)),'g');
grid on;
ylim([1e-7,1])
ylabel('\epsilon')
xlabel('R')
legend('NA', 'RCUs','RCUs SP', 'outage','Is')
end
%-------------------------------------
% SAVE DATA
%-------------------------------------
if DEBUG == 0
    data.M = M;
    data.SNR=snr_db;
    data.n=n;
    data.b=b_vec;
    data.R=R_vec;
    data.eps_out = eps_out;
    data.eps_na = eps_na;
    data.eps_rcus = eps_rcus;
    data.eps_rcus_sp = eps_rcus_sp;
    
    filename = ['R_vs_error_SNR_' num2str(snr_db) '_M_' num2str(M) '.mat'];
    save(filename, 'data', '-v7.3');
    csvwrite(['section2_na_' filename '.csv'], [R_vec' , eps_na']);
    csvwrite(['section2_rcus_' filename '.csv'], [R_vec' , eps_rcus']);
    csvwrite(['section2_rcus_sp_' filename '.csv'], [R_vec' , eps_rcus_sp']);
    csvwrite(['section2_outage_' filename '.csv'], [R_vec' , eps_out']);
    csvwrite(['section2_Is_' filename '.csv'], [ones(length(R_vec),1)*I1 , linspace(1e-7,1,length(R_vec))']);
end
end


function eps_rcus = getErrorProbability(s, n, rho, R,  g_list ,ghat_list,sigma_sq_list)
% Function getErrorProbability(s, n, rho, R,  g_list ,ghat_list,sigma_sq_list)
% that computes the quantities needed to compute the saddlepoint
% approximation and computes the saddlepoint approximation. 
nbrOfRealizations = length(g_list);
eps_ul = nan(1,nbrOfRealizations);
% Implementation of the saddlepoint approximation:
for j = 1:nbrOfRealizations
    g = g_list(j);
    ghat = ghat_list(j);
    sigma_sq = sigma_sq_list(j);
    % Parameters needed to obtain the CGF and then compute saddlepoint:
    betaA_ul = s*rho*abs(g-ghat)^2 + s*sigma_sq;
    betaB_ul = s*(rho*abs(g)^2 + sigma_sq) / (1+s*rho*abs(ghat)^2);
    
    sigma_v = abs(g)^2 *rho + sigma_sq;
    gamma = s/(1 + s*rho*abs(ghat)^2);
    nu_ul = s*gamma*abs(sigma_v - rho* g'*ghat)^2 / (betaA_ul*betaB_ul);
    
    preterm_ul = log(1+s*rho * abs(ghat)^2);
    % Saddlepoint approximation:
    eps_ul(j) = saddlepoint_approximation(n, R, betaA_ul, betaB_ul, nu_ul, preterm_ul);
end
eps_rcus = mean(eps_ul);
end

function H1 = generateChannel(M, nbrOfRealizations)
R1 = eye(M);
H1 = sqrt(0.5)*sqrtm(R1)*(randn(M,nbrOfRealizations)+1i*randn(M,nbrOfRealizations));
end

function [nz, nv, ng] = generateRVs(rho, n, M, nbrOfRealizations)
% Function generateRVs(rho, n, M, nbrOfRealizations) that generates
% the random variables needed in the computation of the generalized
% information density
nz = nan(1,nbrOfRealizations);
nv= nan(1,nbrOfRealizations);
ng= nan(1,nbrOfRealizations);
for i = 1:nbrOfRealizations
    H =  sqrt(0.5)*(randn(M,1)+1i*randn(M,1));
    nH = norm(H);
    g = H'*H / nH;
    
    z = sqrt(0.5)*(randn(M,n)+1i*randn(M,n));
    z = H'*z/nH;
    
    q = sqrt(0.5)*sqrt(rho)*(randn(1,n)+1i*randn(1,n));
    v = g*q + z;
    
    nz(i) = norm(z)^2;
    nv(i)=norm(v)^2;
    ng(i) = abs(g)^2;
end
end



function eps_na = getErrorProbabilityNA(s, n, rho, b,  g_list ,ghat_list,sigma_sq_list)
% Function getErrorProbability(s, n, rho, R,  g_list ,ghat_list,sigma_sq_list)
% that computes the quantities needed to compute the CGFs and its
% derivatives to then compute Is and Vs and compute the normal
% approximation
nbrOfRealizations = length(g_list);
eps_ul = nan(1,nbrOfRealizations);
% Implementation of the saddlepoint approximation:
for j = 1:nbrOfRealizations
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
    % Compute Is and Vs:
    Is = preterm_ul + betaA_ul - betaB_ul;
    Vs = (betaA_ul - betaB_ul).^2 + 2*betaA_ul*betaB_ul*(1-nu_ul);
    % Normal approximation:
    eps_ul(j) = qfunc((n*Is-log(2^b-1)) ./ sqrt(n*Vs));
end
eps_na = mean(eps_ul);
end