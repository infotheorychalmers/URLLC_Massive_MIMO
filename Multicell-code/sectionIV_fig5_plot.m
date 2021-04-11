% Script to process the saved files if the function sectionIV_fig5.m has
% been used.
% This script needs manipulation to process the type of data that wants to
% be generated
clear
%FLAGS (Adapt according to the desired parameters to plot)
np = [10,20,40,60,80,100,120,140,160,180,200,230,250]; %number of pilots
n = 300; %total number of channle uses
b = 8*20; %transmission rate
snr_db = 10; %transmit power [dBm] (10,15,20)
Mlist = [100]; %the number of antennas considered
COMBINER = 'MR'; %what combiner to use [MR, M-MMSE,RZF]
ESTIMATOR = 'MMSE'; %what combiner to use [LS, MMSE]
antennaType = 'ULA';

ASD=25; % angular delay spread
K=5; % number of UEs per cell
L=4; % number of cells

eps_target = 1e-5;

figure
ul_availability = nan(length(Mlist), length(np));
dl_availability = nan(length(Mlist), length(np));

for Mindx = 1:length(Mlist)
    M=Mlist(Mindx);
    for np_indx = 1:length(np)
        filename = [ESTIMATOR '_' COMBINER '_M_' num2str(M) '_SNR_' num2str(snr_db) '_np_' num2str(np(np_indx)) '_n_' num2str(n)  '_ASD_' num2str(ASD) '_K_' num2str(K) '_L_' num2str(L) '.mat'];
        folder = [pwd '/data/'];
        
        simNr = 1:4;
        
        for i = 1:length(simNr)
            cur_file = ['simNr_' num2str(simNr(i)) '_Multicell_' antennaType '_CDF_' filename];
            
            try
                load([folder, cur_file])
            catch
                disp(['The file ' cur_file ' was not found so we skip it.'])
                continue
            end
            
            N = length(data.avg_error_ul);
            
            eps_ul((i-1)*N+1 : i*N) = data.avg_error_ul;
            eps_dl((i-1)*N+1 : i*N) = data.avg_error_dl;
        end
        
        %Obtain UL params
        [cdf_ul, x_ul] = ecdf(eps_ul);
        
        loglog(x_ul, cdf_ul); hold on;
        
        [~ , target_indx] = min(abs(x_ul-eps_target));
        
        if max(x_ul) < eps_target
            target_indx = length(cdf_ul);
        end
        
        ul_availability(Mindx, np_indx) = cdf_ul(target_indx);
        
        
        %Obtain DL params
        [cdf_dl, x_dl] = ecdf(eps_dl);
        
        [~ , target_indx] = min(abs(x_dl-eps_target));
        
        if max(x_dl) < eps_target
            target_indx = length(cdf_dl);
        end
        dl_availability(Mindx, np_indx) = cdf_dl(target_indx);
        
%         subplot(1,2,1)
%         loglog(x_ul,cdf_ul); grid on; xlim([1e-5,1]); hold on; title(['UL, ASD=' num2str(ASD) ', ' antennaType]);
%         
%         subplot(1,2,2)
%         loglog(x_dl,cdf_dl); grid on; xlim([1e-5,1]); hold on; title(['DL, ASD=' num2str(ASD) ', ' antennaType]);
%         
%         
        
%         savename_ul = ['Sec4_UL_' filename '.csv']
%         savename_dl = ['Sec4_DL_' filename '.csv']
%         
%         
%         csvwrite(savename_ul, [x_ul , cdf_ul]);
%         csvwrite(savename_dl, [x_dl , cdf_dl]);
    end
    
end
%legend('np=10','np=20', 'np=40', 'np=60','np=100', 'np=200', 'np=250')
%legend('M=10','M=25', 'M=50', 'M=75','M=100', 'M=125', 'M=150', 'M=175', 'M=200')


% subplot(1,2,1)
% semilogx(Mlist, ul_availability);hold on;
% grid on;
% subplot(1,2,2)
% semilogx(Mlist, dl_availability); hold on;
% grid on;
   
subplot(1,2,1)
plot(np, ul_availability);hold on;
grid on;
subplot(1,2,2)
plot(np, dl_availability); hold on;
grid on;

% savename_ul = ['Multicell_M_UL_Availability_' filename '.csv']
% savename_dl = ['Multicell_M_DL_Availability_' filename '.csv']
% csvwrite(savename_ul, [Mlist' , ul_availability]);
% csvwrite(savename_dl, [Mlist' , dl_availability]);

savename_ul = ['Multicell_np_UL_Availability_' filename '.csv']
savename_dl = ['Multicell_np_DL_Availability_' filename '.csv']
csvwrite(savename_ul, [np' , ul_availability']);
csvwrite(savename_dl, [np' , dl_availability']);
