function section_II_fig1(n,b,snr_db,Mlist,nbrOfRealizations)
% Function section_II_fig1(n,b,snr_db,Mlist,nbrOfRealizations): Generates 
% Fig. 1(a) and Fig. 1(b) in the paper.
% 
% INPUTS:
% n = blocklength
% b = information bits
% snr_db = SNR [dB]
% Mlist = list with number of antennas
% nbrOfRealizations = number of realization in the Monte-Carlo simulations
%
DEBUG = 1;
SNR_normalized = 1;
if DEBUG == 1
    if SNR_normalized == 0
        snr_db = -24;
        Mlist = 50:25:400;
        %Mlist = 50:25:300; % only to compute exact RCUs
    else
        snr_db = 1;
        Mlist = [1 5 10 20 50 100 200 500 1000];
        %Mlist = [1 5 10 20 50 100]; % only to compute exact RCUs
    end
    nbrOfRealizations = 1e3; % Increase to compute the exact RCUs
    n = 100;
    b = 3*20;
end
R = b/n;
% Definition of functions:
C = @(p) log2(1+p); %capacity Gaussian channel
VG = @(p) log2(exp(1))^2 * 2*p./(p+1); %Dispersion Gaussian channel Gaussian inputs
vecnorm = @(A)  sqrt(sum(A.*conj(A),1)); %norm of each column in matrix
% Initializations:
eps_out = nan(1,length(Mlist));
eps_rcus_sp = nan(1,length(Mlist));
eps_rcus = nan(1,length(Mlist));
eps_na = nan(1,length(Mlist));
eps_naG = nan(1,length(Mlist));
% For golden search:
START_INT = 1e-2;
END_INT = 1;
TOL = 1e-2;
for i = 1:length(Mlist)
    disp(['M = ' num2str(Mlist(i))])
    M = Mlist(i);
    rho = 10.^(snr_db/10);%rho = avergage gain * transmit power / sigma^2
    if SNR_normalized == 1
        rho=rho/M; %Normalize by M. Easier visualization
    end
    %-------------------------------------
    % GENERATE CHANNEL REALIZATIONS
    H1 = generateChannel(M, nbrOfRealizations);
    hhat1 = H1;
    v1 = hhat1 ./ (vecnorm(hhat1)); %MR combiner for UE_1: hhat1 = H1;
    sigma_sq = nan(1,nbrOfRealizations);
    g = nan(1,nbrOfRealizations);
    ghat = nan(1,nbrOfRealizations);
    for j=1:nbrOfRealizations
        sigma_sq(j) = v1(:,j)'*v1(:,j); % effective noise variance
        g(j) = v1(:,j)' * H1(:,j) ; % effective channel
        ghat(j) = v1(:,j)' * hhat1(:,j); % effective channel estimate (Perf. CSI)
    end
    %-------------------------------------
    % EVALUATE RCUs for perfect CSI
    f = @(s) getErrorProbability(s, n, rho, R, g,ghat,sigma_sq);     
    [eps_rcus_sp(i), s_val] = golden_search(f, START_INT, END_INT, TOL);
    
    if eps_rcus_sp(i) > 1e-3 % We only compute real RCUs for large epsilon
        [nz,nv,ng]=generateRVs(rho, n, M, nbrOfRealizations); % RVs involved in info. density
        i_s = @(s) -s*nz + (s./(1+s*rho*ng)).*nv + n*log(1+s*rho*ng); % generalized info. density
        error_val = @(s) mean(exp(-max(0, i_s(s) - log(2^b-1)))); % Exact RCUs
        [eps_rcus(i), s_val(i)]=golden_search(error_val,  1e-2, 2, 1e-4); % Optimization over s
    end
    %-------------------------------------
    %EVALUATE OUTAGE
    eps_out(i) = gamcdf((2^R - 1)/rho, M, 1); %compute outage probability
    %-------------------------------------
    %EVALUATE NORMAL APPROXIMATION
    fNA = @(s) getErrorProbabilityNA(s, n, rho, b, g,ghat,sigma_sq); 
    [eps_na(i), s_val_na(i)]=golden_search(fNA,  START_INT, END_INT, TOL); % Optimization over s
    %-------------------------------------
    %EVALUATE NORMAL APPROXIMATION (M->infinity)
    rho_na = M*rho; % Renormalization of M since this is a valid approx. when M->infinity.
    eps_naG(i) =  qfunc((C(rho_na)-R) ./ sqrt(VG(rho_na)./n));
end
%-------------------------------------
% FIGURES
%-------------------------------------
if DEBUG == 1
figure
if SNR_normalized == 0 
    semilogy(Mlist, eps_out); hold on;
    semilogy(Mlist, eps_na);
    semilogy(Mlist, eps_naG);
    semilogy(Mlist, eps_rcus,'black*');
    semilogy(Mlist, eps_rcus_sp,'black');    
else
    loglog(Mlist, eps_out); hold on;
    loglog(Mlist, eps_na);
    loglog(Mlist, eps_naG);
    loglog(Mlist, eps_rcus,'black*');
    loglog(Mlist, eps_rcus_sp,'black');    
end
grid on;
ylim([1e-7,1])
ylabel('\epsilon')
xlabel('M')
legend('Outage','NA','NA ($$M\to\infty$$)', 'RCUs','RCUs SP','Interpreter','latex')
end
%-------------------------------------
% SAVE DATA
%-------------------------------------
if DEBUG == 0
    data.Mlist = Mlist;
    data.SNR=snr_db;
    data.n=n;
    data.b=b;
    data.R=R;
    data.eps_out = eps_out;
    data.eps_na = eps_na;
    data.eps_rcus = eps_rcus;
    data.eps_rcus_sp = eps_rcus_sp;
    filename = ['M_vs_error_SNR_' num2str(snr_db) '.mat'];
    save(filename, 'data', '-v7.3');
    csvwrite(['section2_outage_' filename '.csv'], [Mlist' , eps_out']);
    csvwrite(['section2_na_' filename '.csv'], [Mlist' , eps_na']);
    csvwrite(['section2_rcus_' filename '.csv'], [Mlist' , eps_rcus']);
    csvwrite(['section2_rcus_sp_' filename '.csv'], [Mlist' , eps_rcus_sp']);
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

