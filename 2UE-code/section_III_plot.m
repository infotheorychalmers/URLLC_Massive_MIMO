function section_III_plot


COMBINER = 'MR';
NO_OF_UEs = 2;
angle2 = 40;
np = 2;
n = 300;
ASD=25;


PILOT_CONTAMINATION = 1;



%figure('Units', 'centimeters','Position', [10, 10, 30, 20])



UNCORRELATED = 0;

% BS_CSI = 1;UE_CSI = 1;
% name = createFilename(FIG_NR, COMBINER, NO_OF_UEs, angle2, np, n,ASD, UNCORRELATED, PILOT_CONTAMINATION, BS_CSI, UE_CSI);
% load(name);
% loglog(data.Mlist, data.avg_error_dl, 'r'); hold on;

% BS_CSI = 0;UE_CSI = 2;
% name = createFilename(FIG_NR, COMBINER, NO_OF_UEs, angle2, np, n,ASD, UNCORRELATED, PILOT_CONTAMINATION, BS_CSI, UE_CSI);
% load(name);
% loglog(data.Mlist, data.avg_error_dl, 'b'); hold on;
% 
% BS_CSI = 1;UE_CSI = 0;
% name = createFilename(FIG_NR, COMBINER, NO_OF_UEs, angle2, np, n,ASD, UNCORRELATED, PILOT_CONTAMINATION, BS_CSI, UE_CSI);
% load(name);
% loglog(data.Mlist, data.avg_error_dl, 'mag'); hold on;

BS_CSI = 0;UE_CSI = 0;
name = createFilename(COMBINER, NO_OF_UEs, angle2, np, n,ASD,UNCORRELATED, PILOT_CONTAMINATION, BS_CSI, UE_CSI);
load(name);

subplot(1,2,1)
loglog(data.Mlist, data.avg_error_ul); hold on;
grid on;
ylim([1e-5,1]);

subplot(1,2,2)
loglog(data.Mlist, data.avg_error_dl); hold on;
grid on;
ylim([1e-5,1]);

xlabel('M');

title(['np=' num2str(np)])
legend('DL, UE CSI', 'DL, UE est CSI', 'DL, BS CSI, hardening', 'DL, BS est CSI, hardening')

filename = [COMBINER];
if PILOT_CONTAMINATION
   filename = [filename '_PilotCont'];
end
filename=[filename '.csv'];
csvwrite(['Sec3_M_UL_' filename], [data.Mlist', data.avg_error_ul])
csvwrite(['Sec3_M_DL_' filename], [data.Mlist', data.avg_error_dl])
end


function name = createFilename(COMBINER, NO_OF_UEs, angle2, np, n, ASD, UNCORRELATED, PILOT_CONTAMINATION, BS_CSI, UE_CSI)

folder = [pwd '/Simulation Data/Sec 3/Final/'];

filename = [COMBINER '_UEs_' num2str(NO_OF_UEs) '_angle2_' num2str(angle2) '_np_' num2str(np) '_n_' num2str(n) '_ASD_' num2str(ASD) '.mat'];

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
name = [folder, 'ULA_' filename];

end