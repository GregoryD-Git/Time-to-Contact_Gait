%% Analytical Step Length Difference Calculation - Finley et al. (2015) Neurorehabilitation and Neural Repair
function [SLdiffA,tcomps,scomps,vcomps] = calc_FallRisk_SLdiff_Analytical(matfile,ii,plt,p)

load(matfile)
S = matfile;
fq = eval([S '.FrameRate']);

if p == 'R'
    Vprt = dataKIN.AvgVr;       Vuprt = dataKIN.AvgVl;
    PST = dataKIN.RST;          UST = dataKIN.LST;
    Pffp = dataKIN.Rffp;        Uffp = dataKIN.Lffp;
else
    Vprt = dataKIN.AvgVl;       Vuprt = dataKIN.AvgVr;
    PST = dataKIN.LST;          UST = dataKIN.RST;
    Pffp = dataKIN.Lffp;        Uffp = dataKIN.Rffp;
end

msiz = min([length(Vprt),length(Vuprt),length(PST),length(UST),...
    length(Pffp),length(Uffp)]);
% Inputs require both belt velocities, both step times, and both forward
% foot placements, each as a vector of the same length

%% Calculating step length difference
tcomps = ((Vprt(1:msiz) + Vuprt(1:msiz))./2).*(PST(1:msiz) - UST(1:msiz));
scomps = Uffp(1:msiz) - Pffp(1:msiz);
vcomps = ((PST(1:msiz) + UST(1:msiz))/2).*(Vprt(1:msiz) - Vuprt(1:msiz));

SLdiffA = tcomps + scomps + vcomps;

siz = min([length(SLdiffA),length(dataKIN.LSL),length(dataKIN.RSL)]);

SLdiffi = (dataKIN.LSL(1:siz) - dataKIN.RSL(1:siz))./0.5.*(dataKIN.LSL(1:siz) + dataKIN.RSL(1:siz));

if plt == 1
    figure(ii); hold on
%     subplot(211); hold on
    plot(tcomps,'b')
    plot(scomps,'r')
    plot(vcomps,'m')
    plot(SLdiffA,'k','LineWidth',2)
    ylabel('Analytical SL diff')
%     subplot(212); hold on
%     plot(SLdiffA,SLdiffi,'or')
end

save(matfile,'dataKIN','-append')