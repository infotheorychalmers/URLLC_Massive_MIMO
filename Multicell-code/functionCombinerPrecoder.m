function [V,W] = functionCombinerPrecoder(Hhat,C,R,nbrOfRealizations,M,K,L,p, combiner)
%Compute UL SE for different receive combining schemes using Theorem 4.1.
%
%INPUT:
%Hhat              = M x nbrOfRealizations x K x L x L matrix with the MMSE
%                    channel estimates
%C                 = M x M x K x L x L matrix with estimation error
%                    correlation matrices when using MMSE estimation
%R                 = M x M x K x L x L matrix with spatial correlation
%                    matrices
%tau_c             = Length of coherence block
%tau_p             = Length of pilot sequences
%nbrOfRealizations = Number of channel realizations
%M                 = Number of antennas per BS
%K                 = Number of UEs per cell
%L                 = Number of BSs and cells
%p                 = Uplink transmit power per UE (same for everyone)
%
%OUTPUT:
%SE_MR    = K x L matrix where element (k,l) is the uplink SE of UE k in
%           cell l achieved with MR combining
%SE_RZF   = Same as SE_MR but with RZF combining
%SE_MMMSE = Same as SE_MR but with M-MMSE combining
%SE_ZF    = Same as SE_MR but with ZF combining
%SE_SMMSE = Same as SE_MR but with S-MMSE combining
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
vecnorm = @(A)  sqrt(sum(A.*conj(A),1)); %norm of each column in matrix


%Store identity matrices of different sizes
eyeK = eye(K);
eyeM = eye(M);


if strcmp(combiner, 'M-MMSE')
    %Compute sum of all estimation error correlation matrices at every BS
    C_totM = reshape(p*sum(sum(C,3),4),[M M L]);
end


if strcmp(combiner, 'S-MMSE')
    %Compute sum of intra-cell estimation error correlation matrices at every BS
    CR_totS = zeros(M,M,L);
    for j = 1:L
        CR_totS(:,:,j) = p*(sum(C(:,:,:,j,j),3)+sum(sum(R(:,:,:,[1:j-1 j+1:end],j),3),4));
    end
end



%% Go through all channel realizations
V = nan(M, nbrOfRealizations, K,L);
W = V;
for n = 1:nbrOfRealizations
    
    %Go through all cells
    for j = 1:L
        
        %Extract channel estimate realizations from all UEs to BS j
        Hhatallj = reshape(Hhat(:,n,:,:,j),[M K*L]);
        
        %Compute MR combining in (4.11)
        V_MR = Hhatallj(:,K*(j-1)+1:K*j);
        V(:,n,:,j) = V_MR;
        
        if strcmp(combiner, 'RZF') %Compute RZF combining in (4.9)
            V(:,n,:,j) = p*V_MR/(p*(V_MR'*V_MR)+eyeK);
        end
        
        if strcmp(combiner, 'M-MMSE') %Compute M-MMSE combining in (4.7)
            V(:,n,:,j) = p*((p*(Hhatallj*Hhatallj')+C_totM(:,:,j)+eyeM) \ V_MR); % A \ B = inv(A)*B,   A/B = A*inv(B)
        end
        
        
        if strcmp(combiner, 'ZF') %Compute ZF combining in (4.10)
            V(:,n,:,j) = V_MR/(V_MR'*V_MR+1e-12*eyeK);
        end
        
        W(:,n,:,j) = V(:,n,:,j) ./ vecnorm(V(:,n,:,j));
    end

end
