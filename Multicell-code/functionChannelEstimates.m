function [Hhat,C,R,H] = functionChannelEstimates(R,channelGaindB,nbrOfRealizations,M,K,L,p,tau_p,ESTIMATOR)
%Generate the channel realizations and estimates of these channels for all
%UEs in the entire network. The channels are assumed to be correlated
%Rayleigh fading. The MMSE estimator, EW-MMSE estimator, and LS estimator
%are used. The latter two estimators are only computed if their outputs are
%requested when calling the function.
%
%INPUT:
%R                 = M x M x K x L x L matrix with spatial correlation
%                    matrices for all UEs in the network. R(:,:,k,j,l) is
%                    the correlation matrix for the channel between UE k
%                    in cell j and the BS in cell l. This such matrix can
%                    either include the average channel gain or can be
%                    normalized arbitrarily.
%channelGaindB     = K x L x L matrix containing the average channel gains
%                    in dB of all the channels, if these are not already
%                    included in the spatial channel correlation matrices.
%                    The product R(:,:,k,j,l)*10^(channelGaindB(k,j,l)/10)
%                    is the full spatial channel correlation matrix.
%nbrOfRealizations = Number of channel realizations
%M                 = Number of antennas per BS
%K                 = Number of UEs per cell
%L                 = Number of BSs and cells
%p                 = Uplink transmit power per UE (same for everyone)
%f                 = Pilot reuse factor
%ESIMATOR          = channel estimation scheme
%OUTPUT:
%Hhat_MMSE    = M x nbrOfRealizations x K x L x L matrix with the MMSE
%               channel estimates. The matrix Hhat_MMSE(:,n,k,j,l) is the
%               n:th channel estimate of the channel between UE k in cell j
%               and the BS in cell l.
%C_MMSE       = M x M x K x L x L matrix with estimation error correlation
%               matrices when using MMSE estimation. The matrix is
%               organized in the same way as R.
%tau_p        = Length of pilot sequences
%R            = Scaled version of the input spatial correlation matrices R,
%               where the channel gains from channelGaindB are included
%H            = M x nbrOfRealizations x K x L x L matrix with the true
%               channel realizations. The matrix is organized in the same
%               way as Hhat_MMSE.
%Hhat_EW_MMSE = Same as Hhat_MMSE, but using the EW-MMSE estimator
%C_EW_MMSE    = Same as C_MMSE, but using the EW-MMSE estimator
%Hhat_LS      = Same as Hhat_MMSE, but using the LS estimator
%C_LS         = Same as C_MMSE, but using the LS estimator
%
%
%This Matlab function was developed to generate simulation results to:
%
%Emil Bjornson, Jakob Hoydis and Luca Sanguinetti (2017),
%"Massive MIMO Networks: Spectral, Energy, and Hardware Efficiency",
%Foundations and Trends in Signal Processing: Vol. 11, No. 3-4,
%pp. 154-655. DOI: 10.1561/2000000093.
%
%For further information, visit: https://www.massivemimobook.com
%
%This is version 1.0 (Last edited: 2017-11-04)
%
%License: This code is licensed under the GPLv2 license. If you in any way
%use this code for research that results in publications, please cite our
%monograph as described above.


%% Generate channel realizations

%Generate uncorrelated Rayleigh fading channel realizations
H = (randn(M,nbrOfRealizations,K,L,L)+1i*randn(M,nbrOfRealizations,K,L,L));

%Prepare a matrix to save the channel gains per UE
betas = zeros(K,L,L);


%Go through all channels and apply the channel gains to the spatial
%correlation matrices
for j = 1:L
    
    for l = 1:L
        
        for k = 1:K
            
            if channelGaindB(k,j,l)>-Inf
                
                %Extract channel gain in linear scale
                betas(k,j,l) = 10^(channelGaindB(k,j,l)/10);
                
                %Apply channel gain to correlation matrix
                R(:,:,k,j,l) = betas(k,j,l) * R(:,:,k,j,l);
                
                %Apply correlation to the uncorrelated channel realizations
                Rsqrt = sqrtm(R(:,:,k,j,l));
                H(:,:,k,j,l) = sqrt(0.5)*Rsqrt*H(:,:,k,j,l);
                
            else
                
                betas(k,j,l) = 0;
                R(:,:,k,j,l) = 0;
                H(:,:,k,j,l) = 0;
                
            end
            
        end
        
    end
    
end



%% Perform channel estimation

%setdiff
%union
%assign pilots to the UEs
% This works only for L=1 and L=4!
pilotPattern = nan(K,L);

for i = 1:L
    if i == 1
        pilotSequences = 1:tau_p;
    elseif i ==2 
        pilotSequences = setdiff(1:tau_p, pilotPattern(:,1)); %begin looking for sequences in P \ P1
    elseif i == 3
        pilotSequences = setdiff(1:tau_p, union(pilotPattern(:,1), pilotPattern(:,2))); %begin looking for sequences in P\(P1 U P2)
    else
        pilotSequences = setdiff(1:tau_p, union(union(pilotPattern(:,1), pilotPattern(:,2)) ,pilotPattern(:,3))); %begin looking for sequences in P\(P1 U P2 U P3)
    end
    for j = 1:K
        if isempty(pilotSequences)
            if i == 2 %try to find sequences in P \ P2
                pilotSequences = setdiff(1:tau_p, pilotPattern(:,2)); %create new available patterns such that there will be no duplicate in the cell
            elseif i == 3
                    pilotSequences = setdiff(1:tau_p, union(pilotPattern(:,1),pilotPattern(:,3)));% try to find sequences in P \ (P3 U P1)
                    if isempty(pilotSequences) %if still no sequences
                        pilotSequences = setdiff(1:tau_p, union(pilotPattern(:,2), pilotPattern(:,3)));% try to find sequences in P \ (P3 U P2)
                        if isempty(pilotSequences) %if still no sequences
                            pilotSequences = setdiff(1:tau_p,  pilotPattern(:,3));% try to find sequences in P \ P3
                        end
                    end
            elseif i == 4
                pilotSequences = setdiff(1:tau_p, union(union(pilotPattern(:,2),pilotPattern(:,3)), pilotPattern(:,4)));% try to find sequences in P \ (P2 U P3 U P4)
                if isempty(pilotSequences) %if still no pilot sequences
                    pilotSequences = setdiff(1:tau_p, union(pilotPattern(:,2), pilotPattern(:,4)));% try to find sequences in P \ (P2 U P4)
                    if isempty(pilotSequences)
                        pilotSequences = setdiff(1:tau_p, union(pilotPattern(:,3), pilotPattern(:,4)));% try to find sequences in P \ (P3 U P4)
                        if isempty(pilotSequences)
                            pilotSequences = setdiff(1:tau_p,  union(pilotPattern(:,1), pilotPattern(:,4)));% try to find sequences in P \  (P1 U P4)
                            if isempty(pilotSequences)
                                pilotSequences = setdiff(1:tau_p,   pilotPattern(:,4));% try to find sequences in P \ P4
                            end
                        end
                    end
                end
            end
        end
    
        seqNr = datasample(pilotSequences,1,'Replace',false); %sample a pilot sequence
        pilotPattern(j,i) = seqNr; %assign UE j the sequence
        pilotSequences = pilotSequences(pilotSequences ~= seqNr); %remove pilot sequence from available sequences
    end
end
% [uv,~,idx] = unique(pilotPattern);
% pilotNumbers = accumarray(idx(:),1) %see how many of each sequence we have
% a = sum(pilotNumbers == 1:4)
%%Store identity matrix of size M x M
eyeM = eye(M);

%Prepare for estimation

%Prepare to store channel estimates
Hhat = zeros(M,nbrOfRealizations,K,L,L);

%Prepare to store estimation error correlation matrices
C = zeros(M,M,K,L,L);



%% Go through all cells
for j = 1:L
    Np = sqrt(0.5)*(randn(M,nbrOfRealizations) + 1i*randn(M,nbrOfRealizations)); %same Np for all receptions at the BS
    for l = 1:L
        %Go through all UEs
        for k = 1:K
            
            groupMembers = pilotPattern == pilotPattern(k,l); %Find all UEs that share the same pulot sequence as UE k in cell l
            [ue_share, cell_share] = find(groupMembers == 1); %find out what UEs in what cells use the same pilots
            P = [ue_share, cell_share]; %UEs in same pilot group
            nbrOfSharingUEs = size(P,1);
            
            %Compute processed pilot signal for all UEs that use these pilots, according to (3.5)
            yp_jli=0;
            PsiInv = zeros(M,M);
            for g = 1:nbrOfSharingUEs %go through all the UEs that share the same pilot sequence as UE k in cell l
                UEindx = P(g,1);
                cellindx= P(g,2);
                yp_jli = yp_jli + sqrt(p)*tau_p*H(:,:,UEindx,cellindx,j); %sum all the channels from the UEs tat share the sequence
                PsiInv = PsiInv + p*tau_p*R(:,:,UEindx,cellindx,j) ;%Compute the matrix that is inverted in the MMSE estimator
            end
            %Generate realizations of normalized noise
            yp_jli = yp_jli +  sqrt(tau_p)*Np;
            PsiInv = PsiInv + eyeM;
            
            %Compute a vector with elements that are inverted in the EW-MMSE estimator
            PsiInvDiag = diag(PsiInv);
                        
            if strcmp(ESTIMATOR, 'MMSE')
                %Compute MMSE estimate of channel between BS l and UE k in
                %cell j using (3.9) in Theorem 3.1
                RPsi = R(:,:,k,l,j) / PsiInv;
                Hhat(:,:,k,l,j) = sqrt(p)*RPsi*yp_jli;
                
                %Compute corresponding estimation error correlation matrix, using (3.11)
                C(:,:,k,l,j) = R(:,:,k,l,j) - p*tau_p*RPsi*R(:,:,k,l,j);
            elseif strcmp(ESTIMATOR, 'EW-MMSE')
                %Compute EW-MMSE estimate of channel between BS l and
                %UE k in cell j using (3.33)
                A_EW_MMSE = diag(sqrt(p)*diag(R(:,:,k,l,j)) ./ PsiInvDiag);
                Hhat(:,:,k,l,j) = A_EW_MMSE*yp_jli;
                
                %Compute corresponding estimation error correlation
                %matrix, using the principle from (3.29)
                productAR = A_EW_MMSE * R(:,:,k,l,j);
                
                C(:,:,k,l,j) = R(:,:,k,l,j) - (productAR + productAR') * sqrt(p)*tau_p + tau_p*A_EW_MMSE*PsiInv*A_EW_MMSE';
                
            elseif strcmp(ESTIMATOR, 'LS')
                %Compute LS estimate of channel between BS l and UE k
                %in cell j using (3.35) and (3.36)
                A_LS = 1/(sqrt(p)*tau_p);
                Hhat(:,:,k,l,j) = A_LS*yp_jli;
                
                %Compute corresponding estimation error correlation
                %matrix, using the principle from (3.29)
                productAR = A_LS * R(:,:,k,l,j);
                
                C(:,:,k,l,j) = R(:,:,k,l,j) - (productAR + productAR') * sqrt(p)*tau_p + tau_p*A_LS*PsiInv*A_LS';
                
            else
                Hhat(:,:,k,l,j) = nan;
                C(:,:,k,l,j) = nan;
            end
        end
    end
end
