% Script to process the saved files if the function sectionIII_fig3.m has
% been used. 

%FLAGS (Adapt according to the desired parameters to plot)
np = 2; %number of pilots
n = 300; %total number of channle uses
b = 8*20; %information bits
rho_db = 10; %transmit power [dBm]
Mlist = [100]; %the number of antennas considered
COMBINER = 'M-MMSE'; %what combiner to use [MR, M-MMSE]
PILOT_CONTAMINATION = 0; %Let the two UEs use the same pilot sequence
NO_OF_UEs = 2; % 1 or 2
UNCORRELATED = 0; %0 = spatial correlation
BS_CSI = 0;
UE_CSI = 0; %0 = UE uses E[H1'W1], 1 = UE knows H1'W1, 2 = UE knows W1 and Hhat1,
antennaType = 'ULA';

% Recover file name depending on the parameters:
filename = [COMBINER '_M_' num2str(Mlist) '_UEs_' num2str(NO_OF_UEs) '_SNR_' num2str(rho_db) '_np_' num2str(np) '_n_' num2str(n) '.mat'];

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
filename = [ '_' antennaType '_CDF_' filename];
folder = [pwd '/']; % adapt this according to the folder where the data is saved

eps_ul = nan(1e3,1);
eps_dl = nan(1e3,1);

% Concatenate the error prob. values of all the saved simulations:
simNr = 4;
for i = 1:simNr
    cur_file = ['simNr_' num2str(i), filename];
    try
    load([folder, cur_file])
    catch
       continue 
    end
    N = length(data.avg_error_ul);
    
    eps_ul((i-1)*N+1 : i*N) = data.avg_error_ul;
    eps_dl((i-1)*N+1 : i*N) = data.avg_error_dl;
end

% Obtain the cdfs according to the obtained error probability values
[cdf_ul, x_ul] = ecdf(eps_ul);
x_ul(end+1)=1;
cdf_ul(end+1)=1;

[cdf_dl, x_dl] = ecdf(eps_dl);
x_dl(end+1)=1;
cdf_dl(end+1)=1;

subplot(1,2,1)
loglog(x_ul, cdf_ul); hold on
grid on
xlim([1e-5,1])

subplot(1,2,2)
loglog(x_dl, cdf_dl); hold on
grid on
xlim([1e-5,1])

% Save the CDFs in csv mode:
if PILOT_CONTAMINATION == 1
    savename_ul = ['Sec3_UL_CDF_iCSI' num2str(COMBINER) '_pilot_cont_M_' num2str(Mlist) '_' antennaType '.csv']
    savename_dl = ['Sec3_DL_CDF_iCSI' num2str(COMBINER) '_pilot_cont_M_' num2str(Mlist) '_' antennaType '.csv']
else
    savename_ul = ['Sec3_UL_CDF_iCSI' num2str(COMBINER) '_M_' num2str(Mlist) '_' antennaType '.csv']    
    savename_dl = ['Sec3_DL_CDF_iCSI' num2str(COMBINER) '_M_' num2str(Mlist) '_' antennaType '.csv']
end

csvwrite(savename_ul, [x_ul , cdf_ul]);
csvwrite(savename_dl, [x_dl , cdf_dl]);





