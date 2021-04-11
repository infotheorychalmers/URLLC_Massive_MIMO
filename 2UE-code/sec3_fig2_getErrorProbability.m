function [avg_error, s_val, epsilon] = sec3_fig2_getErrorProbability(n, rho, b, g_list, g_hat_list, sigma_sq_list, nbrOfRealizations, s_start)
% Function sec3_getErrorProbability(n, rho, b, g_list ,ghat_list,sigma_sq_list, nbrOfRealizations, s_start))
% that computes the saddlepoiint approximation of the error probability given by the RCUs
if nargin < 8
   s_start = 1; 
end

rate = b/n; % rate in bits

% Function that computes the saddlepoint as a function of the parameter s
f_avg_error = @(s) sec3_SaddlepointApprox(s, n, rho, rate, g_list,g_hat_list,sigma_sq_list, nbrOfRealizations);
% Function to optimize over s
[avg_error,  s_val, epsilon] = searchForCandidateS(f_avg_error,1e-10,s_start);

end

function epsilon = sec3_SaddlepointApprox(s, n, rho, rate, g_list, g_hat_list, sigma_sq_list, nbrOfRealizations)
% Function sec3_SaddlepointApprox(s, n, rho, rate, g_list, g_hat_list, sigma_sq_list, nbrOfRealizations)
% that computes the quantities needed to compute the CGFs and its
% derivatives to then compute saddlepoiint approximation of the
% error probability given by the RCUs

epsilon = nan(1, nbrOfRealizations);

for j = 1:nbrOfRealizations
    % Get channel, channel estimate, and effective noise:
    g = g_list(j);
    ghat = g_hat_list(j);
    sigma_sq = sigma_sq_list(j);
    % Parameters related to the CGF of the info density:
    betaA_ul = s*rho*abs(g-ghat)^2 + s*sigma_sq;
    betaB_ul = s*(rho*abs(g)^2 + sigma_sq) / (1+s*rho*abs(ghat)^2);
    sigma_v = abs(g)^2 *rho + sigma_sq;
    gamma = s/(1 + s*rho*abs(ghat)^2);
    nu_ul = s*gamma*abs(sigma_v - rho* g'*ghat)^2 / (betaA_ul*betaB_ul);
    preterm_ul = log(1+s*rho * abs(ghat)^2);
    % Compute the saddlepoint approximation
    epsilon(j) = saddlepoint_approximation(n, rate, betaA_ul, betaB_ul, nu_ul, preterm_ul);
end
end

function [error_prob, s_val, epsilon] = searchForCandidateS(f, eps_target, s_start)
% Function to optimize over the parameter s:

s_list = fliplr(logspace(-8,log10(s_start),50));
s_list(s_list > s_start) = [];
eps_debug = [];

eps_old = inf;
epsilon_old = [];
for ii = 1:length(s_list)
    s_candidate = s_list(ii);
    
    epsilon = f(s_candidate);
    eps_cur = mean(epsilon);
    if eps_cur < eps_target
        break;
    end
    
    eps_debug(ii)=eps_cur;
    
    if eps_cur > eps_old
        eps_cur = eps_old;
        epsilon = epsilon_old;
        s_candidate = s_list(ii-1);
        break;
    else
        eps_old = eps_cur;
        epsilon_old = epsilon;
    end
end

error_prob = eps_cur;
s_val = s_candidate;
end


