%% Time to Contact of swing foot and center of mass to physical and anticipated boundaries
function calc_TtC_Ka_Kp(matfile,pp,pltCHK,pltSav,KaLim)

load(matfile,matfile,'Events','dataKIN','data3D','dataVel','CoMoffset',...
    'Lspeed','Rspeed','Slowlimb')

kfq = eval([matfile '.FrameRate']);             % kinematic frame rate
afq = 1000;      % analog sample rate
scf = 10; % scaling factor between kinematic and analog
n = 101; % interp1 scaled to 0-100% swing phase
%% ---------------------------- Get data
% toe-off and heel strike events
R_Toff = Events.RTOidx(Events.RTOidx<(Events.RICidx(end))); rt1 = R_Toff(1); 
R_Hstk = Events.RICidx(Events.RICidx>rt1);
L_Toff = Events.LTOidx(Events.LTOidx<(Events.LICidx(end))); lt1 = L_Toff(1); 
L_Hstk = Events.LICidx(Events.LICidx > lt1);

R_ToffK = round(R_Toff/scf);     R_HstkK = round(R_Hstk/scf);
L_ToffK = round(L_Toff/scf);     L_HstkK = round(L_Hstk/scf);
LEvent = length(L_Toff);         REvent = length(R_Toff);
EvLen = min([LEvent REvent]);

% ensure first event doesn't equal zero due to rounding erorr
R_ToffK(R_ToffK==0) = 1;        R_HstkK(R_HstkK==0) = 1;
L_ToffK(L_ToffK==0) = 1;        L_HstkK(L_HstkK==0) = 1;

% Center of mass position and velocity
if isfield(data3D,'COM2Y') && ~isnan(data3D.COM2Y(1))
    COMY = data3D.COM2Y+CoMoffset;
else
    COMY = data3D.COM1Y-CoMoffset;
end

dt = data3D.Time(1); t = data3D.Time;
% COM position and velocity, units = meters, meters/sec
CoMp = COMY./1000; cv = diff(CoMp)./dt;    
% Concatenating the last data point to ensure all vectors are the same length
CoMv = [cv;cv(end,:)]; 

% scaling ANKLE/TOE & COP data to meters
LAp = (data3D.LANKLEY)./1000;   RAp = (data3D.RANKLEY)./1000;
LTp = (data3D.LTOEY)./1000;     RTp = (data3D.RTOEY)./1000; %LTp(LTp>0.5) = NaN; RTp(RTp<-0.6 | RTp>0.6) = NaN;

% Ankle velocity - contatenating last data point to preserve length
Lv = diff(LAp)./dt;             Rv = diff(RAp)./dt;
LAv = [Lv;Lv(end)];             RAv = [Rv;Rv(end)];

% if exist('dataCOP','var')
% Lcopap = dataCOP.COPLeftY./1000; Rcopap = dataCOP.COPRightY./1000;
% LCoP = interp1(dataForce.Time,Lcopap,dataForce.Time(1:10:end),'makima');
% RCoP = interp1(dataForce.Time,Rcopap,dataForce.Time(1:10:end),'makima');
% end


%% ---------------------------- Cut data
% Cut data into segments based on toe-off and heel strike timings for each limb
% Each segment represents swing phase from toe-off to heel strike of the
% ipsilateral limb

% preallocating space for non-time normalized data
Cp_pL = NaN(LEvent-1,300);  Cv_pL = Cp_pL;  Sfp_pL = Cp_pL; Stp_pL = Cp_pL;         
Sf_aL = Cp_pL;              SFv_pL = Cp_pL; St_aR = Cp_pL; Bp_pL = Cp_pL; 
Ba_L = Cp_pL;               Bt_pL = Cp_pL;  RA_pL = Cp_pL;

Cp_pR = Cp_pL;              Cv_pR = Cp_pR;  Sfp_pR = Cp_pL; Stp_pR = Cp_pL;         
Sf_aR = Cp_pL;              SFv_pR = Cp_pL; St_aL = Cp_pL; Bp_pR = Cp_pL; 
Ba_R = Cp_pL;               Bt_pR = Cp_pL;  LA_pR = Cp_pL;

% cutting data using toe-off and heel strike events
% left TtC is for left foot in swing, right in stance
% B_p is then the right foot and B_a is the left foot

% shift window to capture xx frames before and xx frames after toe-off and
% heel strike events
sW = [0 0]; % shift Window

%% Left TtC calculations
for ii = 1:EvLen-1
    % ------------ PHYSICAL BOUNDARY CALCS --------------------------------
    % CoMp
    Lcmp = CoMp(L_ToffK(ii)+sW(1):L_ToffK(ii+1)+sW(2)); LcmpLen = length(Lcmp); % L_HstkK(ii)
    Cp_pL (ii,1:LcmpLen) = Lcmp;
    
    % CoMv - must add belt speed to obtain relative speed of CoM
    if exist('dataVel','var')
        LSp = dataVel.LVelf(L_ToffK(ii)+sW(1):L_ToffK(ii+1)+sW(2)); nLVelLen = length(LSp);
        RSp = dataVel.RVelf(L_ToffK(ii)+sW(1):L_ToffK(ii+1)+sW(2)); nRVelLen = length(RSp);
    else
        LSp = Lspeed;
        RSp = Rspeed;
    end
    % CoM adds Bp movement speed to account for relative speed
    Lcmv = CoMv(L_ToffK(ii)+sW(1):L_ToffK(ii+1)+sW(2)); LcmvLen = length(Lcmv);
    Cv_pL(ii,1:LcmvLen) = Lcmv + RSp;
    
    % Ankle - swing foot
    Lank = LAp(L_ToffK(ii)+sW(1):L_ToffK(ii+1)+sW(2)); LApLen = length(Lank);
    Ltoe = LTp(L_ToffK(ii)+sW(1):L_ToffK(ii+1)+sW(2)); LTpLen = length(Ltoe);
    Sfp_pL(ii,1:LApLen) = Lank;
    Stp_pL(ii,1:LTpLen) = Ltoe;
   
    
    % sf adds Bp movement speed to account for relative speed
    LaVel = LAv(L_ToffK(ii)+sW(1):L_ToffK(ii+1)+sW(2)); LAvLen = length(LaVel);
    SFv_pL(ii,1:LAvLen) = LaVel + RSp;

    % CoP - physical boundary
    % Left B_p is the right foot during left swing
    % Left B_a is the left foot at end swing - heel strike
    RA_pB = RAp(L_ToffK(ii)+sW(1):L_ToffK(ii+1)+sW(2)); RApLen = length(RA_pB);
    RT_pB = RTp(L_ToffK(ii)+sW(1):L_ToffK(ii+1)+sW(2)); RTpLen = length(RT_pB);
    RA_pL(ii,1:RApLen) = RA_pB;
    Bt_pL(ii,1:RTpLen) = RT_pB;
    
    % ------------ ANTICIPATED BOUNDARY CALCS -----------------------------
    % anticipated boundary - i.e. next step
    % left foot from left heel strike to next left toe-off
    LA_aB = LAp(L_HstkK(ii)+sW(1):L_ToffK(ii+1)+sW(2)); LAaLen = length(LA_aB);
    LT_aB = LTp(L_HstkK(ii)+sW(1):L_ToffK(ii+1)+sW(2)); LTaLen = length(LT_aB);
%     LA_aB = LAp(L_HstkK(ii)+sW(1):L_HstkK(ii)+sW(2)); LAaLen = length(LA_aB);
%     LT_aB = LTp(L_HstkK(ii)+sW(1):L_HstkK(ii)+sW(2)); LTaLen = length(LT_aB);
    Ba_L(ii,1:LAaLen) = LA_aB;
    % right foot during next left step
    sfaR = RAp(L_HstkK(ii)+sW(1):L_ToffK(ii+1)+sW(2)); SFpaLen = length(sfaR);
    staR = RTp(L_HstkK(ii)+sW(1):L_ToffK(ii+1)+sW(2)); STpaLen = length(sfaR);
%     sfaR = RAp(L_HstkK(ii)+sW(1):L_HstkK(ii)+sW(2)); SFpaLen = length(sfaR);
%     staR = RTp(L_HstkK(ii)+sW(1):L_HstkK(ii)+sW(2)); STpaLen = length(sfaR);
    Sf_aR(ii,1:SFpaLen) = sfaR;
    St_aR(ii,1:STpaLen) = staR;
end

%% Right TtC calculations
for ii = 1:EvLen-1
    % ------------ PHYSICAL BOUNDARY CALCS --------------------------------
    % CoMp
    Rcmp = CoMp(R_ToffK(ii)+sW(1):R_ToffK(ii+1)+sW(2)); RcmpLen = length(Rcmp);
    Cp_pR (ii,1:RcmpLen) = Rcmp;
    
    % CoMv - must add belt speed to obtain relative speed of CoM
    if exist('dataVel','var')
        LSp = dataVel.LVelf(R_ToffK(ii)+sW(1):R_ToffK(ii+1)+sW(2)); nLVelLen = length(LSp);
        RSp = dataVel.RVelf(R_ToffK(ii)+sW(1):R_ToffK(ii+1)+sW(2)); nRVelLen = length(RSp);
    else
        LSp = Lspeed;
        RSp = Rspeed;
    end
    % CoM adds Bp movement speed to account for relative speed
    Rcmv = CoMv(R_ToffK(ii)+sW(1):R_ToffK(ii+1)+sW(2)); RcmvLen = length(Rcmv);
    Cv_pR(ii,1:RcmvLen) = Rcmv + LSp;
    
    % Ankle - swing foot
    Rank = RAp(R_ToffK(ii)+sW(1):R_ToffK(ii+1)+sW(2)); RApLen = length(Rank);
    Rtoe = RTp(R_ToffK(ii)+sW(1):R_ToffK(ii+1)+sW(2)); RTpLen = length(Rtoe);
    Sfp_pR(ii,1:RApLen) = Rank;
    Stp_pR(ii,1:RTpLen) = Rtoe;
   
    % sf adds Bp movement speed to account for relative speed
    RaVel = RAv(R_ToffK(ii)+sW(1):R_ToffK(ii+1)+sW(2)); RAvLen = length(RaVel);
    SFv_pR(ii,1:RAvLen) = RaVel + LSp;

    % CoP - physical boundary
    % Left B_p is the right foot during left swing
    % Left B_a is the left foot at end swing - heel strike
    LA_pB = LAp(R_ToffK(ii)+sW(1):R_ToffK(ii+1)+sW(2)); LApLen = length(LA_pB);
    LT_pB = LTp(R_ToffK(ii)+sW(1):R_ToffK(ii+1)+sW(2)); LTpLen = length(LT_pB);
    LA_pR(ii,1:RApLen) = LA_pB;
    Bt_pR(ii,1:RTpLen) = LT_pB;
    
    % ------------ ANTICIPATED BOUNDARY CALCS -----------------------------
    % anticipated boundary - i.e. next step
    % right ankle from left heel strike to next left toe-off
    RA_aB = RAp(R_HstkK(ii)+sW(1):R_ToffK(ii+1)+sW(2)); RAaLen = length(RA_aB);
    RT_aB = RTp(R_HstkK(ii)+sW(1):R_ToffK(ii+1)+sW(2)); RTaLen = length(RT_aB);
%     RA_aB = RAp(R_HstkK(ii)+sW(1):R_HstkK(ii)+sW(2)); RAaLen = length(RA_aB);
%     RT_aB = RTp(R_HstkK(ii)+sW(1):R_HstkK(ii)+sW(2)); RTaLen = length(RT_aB);
    Ba_R(ii,1:RAaLen) = RA_aB;
    % left ankle during next left step - i.e. Right ankle
    sfaL = LAp(R_HstkK(ii)+sW(1):R_ToffK(ii+1)+sW(2)); SFpaLen = length(sfaL);
    staL = LTp(R_HstkK(ii)+sW(1):R_ToffK(ii+1)+sW(2)); STpaLen = length(sfaL);
%     sfaL = LAp(R_HstkK(ii)+sW(1):R_HstkK(ii)+sW(2)); SFpaLen = length(sfaL);
%     staL = LTp(R_HstkK(ii)+sW(1):R_HstkK(ii)+sW(2)); STpaLen = length(sfaL);
    Sf_aL(ii,1:SFpaLen) = sfaL;
    St_aL(ii,1:STpaLen) = staL;
end

%% Kp, Cm, and Sf relative to Bp
% Bp (stance ankle) adjustment for toe-marker
Stoffset_L = nanmean((Bt_pL - RA_pL),1);        Stoffset_R = nanmean((Bt_pR - LA_pR),1);
Bp_pL(1:size(Bt_pL,1),:) = RA_pL + Stoffset_L;  Bp_pR(1:size(Bt_pR,1),:) = LA_pR + Stoffset_R;

% Kp - CoM and swing foot coupling ratio with physical boundary
SFoffset_L = nanmean(Stp_pL - Sfp_pL,1);        SFoffset_R = nanmean(Stp_pR - Sfp_pR,1);
Sfpoff_pL = Sfp_pL + SFoffset_L;                Sfpoff_pR = Sfp_pR + SFoffset_R;
CmCrossL = Bp_pL - Cp_pL;                       CmCrossR = Bp_pR - Cp_pR;
Tcm_pL = CmCrossL./Cv_pL;                       Tcm_pR = CmCrossR./Cv_pR;
SfCrossL = Bp_pL - Sfpoff_pL;                   SfCrossR = Bp_pR - Sfpoff_pR;
Tsf_pL = SfCrossL./SFv_pL;                      Tsf_pR = SfCrossR./SFv_pR; 
Kp_L = Tsf_pL./Tcm_pL;                          Kp_R = Tsf_pR./Tcm_pR; 
sL = size(Tcm_pL,2);                            sR = size(Tcm_pR,2);
% calc total length of each strides data before cutting spurius data
Kisnan_L = sL - max(cumsum(isnan(Kp_L),2),[],2); Kisnan_R = sR - max(cumsum(isnan(Kp_R),2),[],2); % frame where NaNs start (i.e. end of swing phase

% Cm and Sf %stance crossing Bp
% % finding multiple crosses within a swing period and removing data after
% % first cross by turning to NaN
% [XLI,XLJ] = find(diff(diff(CmCrossL,[],2)>0,[],2)==1); [Lm,LmI] = min(XLI); CmCrossL(Lm,XLJ(LmI):end) = NaN;
% [XRI,XRJ] = find(diff(diff(CmCrossR,[],2)>0,[],2)==1); [Rm,RmI] = min(XRI); CmCrossR(Rm,XRJ(RmI):end) = NaN;

% find Cm/Sf cross data that starts below zero and turn to NaN - since won't be
% calculated anyway and messes up data calcs to get Cm/Sf % swing crossing
% Bp
CmLremove = CmCrossL(:,1)<0; CmCrossL(CmLremove,:) = NaN; % took out find(CmCrossL(:,1)<0), replace if code doesn't work here
CmRremove = CmCrossR(:,1)<0; CmCrossR(CmRremove,:) = NaN;
SfLremove = SfCrossL(:,1)<0; SfCrossL(SfLremove,:) = NaN;
SfRremove = SfCrossR(:,1)<0; SfCrossR(SfRremove,:) = NaN;

% finding where Cm crosses Bp
[CmCLI,CmCLJ] = find(diff((CmCrossL)<0,[],2)==1);   [CmCRI,CmCRJ] = find(diff((CmCrossR)<0,[],2)==1);   
% finding indexes where Cm crosses Bp
[~,CmCLIsort] = sort(CmCLI);                        [~,CmCRIsort] = sort(CmCRI);
% finding frame where Cm data ends
[CmeLI,CmeLJ] = find(diff((CmCrossL)<0,[],2)==-1);  [CmeRI,CmeRJ] = find(diff((CmCrossR)<0,[],2)==-1); 
% fidning frame where Cm data ends
[~,CmeLIsort] = sort(CmeLI);                        [~,CmeRIsort] = sort(CmeRI);
% finding where Cm crosses Bp
[SfCLI,SfCLJ] = find(diff((SfCrossL)<0,[],2)==1);   [SfCRI,SfCRJ] = find(diff((SfCrossR)<0,[],2)==1);   
% finding indexes where Cm crosses Bp
[~,SfCLIsort] = sort(SfCLI);                        [~,SfCRIsort] = sort(SfCRI);
% finding frame where Cm data ends
[SfeLI,SfeLJ] = find(diff((SfCrossL)<0,[],2)==-1);  [SfeRI,SfeRJ] = find(diff((SfCrossR)<0,[],2)==-1); 
% fidning frame where Cm data ends
[~,SfeLIsort] = sort(SfeLI);                        [~,SfeRIsort] = sort(SfeRI);

% --------------- plot this ----------------------
% percent swing where Cm/Sf cross Bp
Cm_crossL = 100*CmCLJ(CmCLIsort)./CmeLJ(CmeLIsort); Cm_crossR = 100*CmCRJ(CmCRIsort)./CmeRJ(CmeRIsort);
Sf_crossL = 100*SfCLJ(SfCLIsort)./SfeLJ(SfeLIsort); Sf_crossR = 100*SfCRJ(SfCRIsort)./SfeRJ(SfeRIsort);
% ------------------------------------------------

% TtC with only positive values 
TcmIL = Tcm_pL<0; Tcm_pL(TcmIL) = NaN;     TcmIR = Tcm_pR<0; Tcm_pR(TcmIR) = NaN;
TsfIL = Tsf_pL<0; Tsf_pL(TsfIL) = NaN;     TsfIR = Tsf_pR<0; Tsf_pR(TsfIR) = NaN;

% Kp, Kp at entry to UE & Kp at toe-off
Kp_L = Tsf_pL./Tcm_pL;                      Kp_R = Tsf_pR./Tcm_pR;  % recalculate Kp without zero crossing TtC
Kframe_L = sL - sum(cumsum(Kp_L<1,2)>0,2);  Kframe_R = sR - sum(cumsum(Kp_R<1,2)>0,2); % frame where Kp passes 1
Kframe_L(Kframe_L==sL) = NaN;               Kframe_R(Kframe_R==sL) = NaN;
% ------------ plot this ----------------
KpXone_L = 100*(Kframe_L./Kisnan_L);        KpXone_R = 100*(Kframe_R./Kisnan_R); % percent swing of Kp less than 1
Kp_tL = Kp_L(:,1);                          Kp_tR = Kp_R(:,1); % Kp at toe-off

Tcm_tL = Tcm_pL(:,1);                       Tcm_tR = Tcm_pR(:,1);
Tsf_tL = Tsf_pL(:,1);                       Tsf_tR = Tsf_pR(:,1);

% Tcm and Tsf values can be missing (NaN) when the cm or sf has already
% passed the physical boundary at the time of trailing limb toe-off,
% therefore NaNs are set to "0" to ensure there is data during missing
% strides
% Tcm_tL(isnan(Tcm_tL)) = 0;                  Tcm_tR(isnan(Tcm_tR)) = 0;
% Tsf_tL(isnan(Tsf_tL)) = 0;                  Tsf_tR(isnan(Tsf_tR)) = 0;

% save data - Kp at toe-off, Tcm at toe-off, Tsf at toe-off
dataKIN.Kp_tL = Kp_tL;                      dataKIN.Kp_tR = Kp_tR;
dataKIN.Tcm_tL = Tcm_tL;                    dataKIN.Tcm_tR = Tcm_tR;
dataKIN.Tsf_tL = Tsf_tL;                    dataKIN.Tsf_tR = Tsf_tR;
% ---------------------------------------

% plot saved parameters - Kp_t, Tcm, Tsf
if pltSav == 1
    figure(pp); clf; 
    subplot(2,3,1); hold on
    p1 = plot(1:length(Kp_tL),Kp_tL,'-k');          plot(1:length(Kp_tR),Kp_tR,'-r');
    legend('left','right');  legend box off;        title('Kp_{Bp}')
    box on;                                         ylabel({'physical boundary','ratio'})
    
    subplot(2,3,2); hold on
    plot(1:length(Tcm_tL),Tcm_tL,'-k');             plot(1:length(Tcm_tR),Tcm_tR,'-r');
    ylabel('time (s)');                             xlabel('strides')
    box on;                                         title('TtC_{cm} at toe-off');    
    
    subplot(2,3,3); hold on
    plot(1:length(Tsf_tL),Tsf_tL,'-k');             plot(1:length(Tsf_tR),Tsf_tR,'-r');
    ylabel('time (s)');                             xlabel('strides')
    box on;                                         title('TtC_{sf} at toe-off');   
    
    sgtitle([matfile(3:7) ' file #:' matfile(11:12) ' slow limb: ' Slowlimb])
end

% plot check all parameters
if pltCHK == 1
    figure(pp); clf; 
    % ------------ plot Kp, Tcm, Tsf ------------------
    % left
    subplot(2,4,1); hold on
    yyaxis left;    p1 = plot(1:size(Kp_L,2),Kp_L,'-r');
    ylabel('Kp');   set(gca,'ycolor','k'); y1 = get(gca,'ylim');
    ylim([0 min([y1(2) 4])]);              x1 = get(gca,'xlim');
    if x1(2) > 50; xlim([0 50]); end
    yyaxis right;   p2 = plot(1:size(Kp_L,2),Tcm_pL,'-k');
    p3 = plot(1:size(Kp_L,2),Tsf_pL,'-b');
    ylabel('TtC');  title('left limb');     ylim([0 1])
    legend([p1(1) p2(1) p3(1)],'Kp','Tcm','Tsf'); legend box off
    box on;         set(gca,'ycolor','k'); y2 = get(gca,'ylim');
    xlabel('frames'); ylim([0 min([y2(2) 4])]); if x1(2) > 50; xlim([0 50]); end
    
    % right
    subplot(2,4,2); hold on
    yyaxis left;    plot(1:size(Kp_R,2),Kp_R,'-r');
    ylabel('Kp');   set(gca,'ycolor','k');  y3 = get(gca,'ylim');
    ylim([0 min([y3(2) 4])]);               x2 = get(gca,'xlim');
    if x2(2) > 50; xlim([0 50]); end
    yyaxis right;   plot(1:size(Kp_R,2),Tcm_pR,'-k');
    plot(1:size(Kp_R,2),Tsf_pR,'-b');           ylim([0 1])
    ylabel('TtC');  title('right limb')
    box on;         set(gca,'ycolor','k');  y4 = get(gca,'ylim');
    xlabel('frames');  ylim([0 y4(2)]);     if x1(2) > 50; xlim([0 50]); end
    
    % ------------ plot stride-by-stride Kp1,Kp_t,Tcm_t, Tsf_t ------------
    subplot(2,4,3); hold on
    yyaxis left
    p4 = plot(1:length(Kp_tL),Kp_tL,'-k');          plot(1:length(Kp_tR),Kp_tR,'--k');
    set(gca,'ycolor','k'); ylabel('Kp @toe-off')
    yyaxis right
    p5 = plot(1:length(KpXone_L),KpXone_L,'-b');    plot(1:length(KpXone_R),KpXone_R,'--b');
    legend([p4(1) p5(1)],'Kp_t','KpXone');          ylim([0 100]); xlabel('strides')
    set(gca,'ycolor','b'); ylabel('%swing at K_{p} = 1'); legend box off
    box on;   title('Kp');                          text(10,90,[{'solid: left'},{['slow limb: ' Slowlimb]}],'FontSize',8)
    
    subplot(2,4,4); hold on
    p6 = plot(1:length(Tcm_tL),Tcm_tL,'-b'); plot(1:length(Tcm_tR),Tcm_tR,'--b');
    p7 = plot(1:length(Tsf_tL),Tsf_tL,'-r'); plot(1:length(Tsf_tR),Tsf_tR,'--r');
    legend([p6 p7],'Tcm','Tsf');    ylabel('time (s)'); xlabel('strides')
    box on;       title('TtC_{Bp} @toe-off'); legend box off
    
end

%% Ka, Cm, and Sf relative to Ba
% anticipated boundary - positions scaled relative to opposite foot
% position to ensure all positions move in the same reference frame

% need to add toe marker offset to ankle marker
% stance foot offset
SFoffset_L = nanmean(St_aR - Sf_aR,1);        SFoffset_R = nanmean(St_aL - Sf_aL,1);
% swing foot offset
SfOff_aR = Sf_aR + SFoffset_L;                SfOff_aL = Sf_aL + SFoffset_R;

% TtC swing foot with Ba
% all positions are relative to the contralateral foot during stance and/or
% swing to account for the Cm never 'reaching' the Ba position
% Ba at left heel strike relative to contralateral swing foot
BaLhs_pL = Ba_L - SfOff_aR;                     BaRhs_pR = Ba_R - SfOff_aL;   % anticipated boundary (i.e heel strike position) relative to contralateral swing foot
% Ipsilateral swing foot relative to Ba moving position relative to Ba
% location at heel strike
SfRel_pL = Sfp_pL - Bp_pL - BaLhs_pL(:,1);      SfRel_pR = Sfp_pR - Bp_pR - BaRhs_pR(:,1);  % swing foot relative to stance foot relative to Ba
SfRel_vL = [diff(SfRel_pL,[],2),diff(SfRel_pL(:,end-1:end),[],2)]./dt;      
SfRel_vR = [diff(SfRel_pR,[],2),diff(SfRel_pR(:,end-1:end),[],2)]./dt;

% % % % % Need to remove position data after the foot passes the Ba
% % % % % The SfRel position 
% % % % % diff(i,[],2)>0
Tsf_aL = abs(SfRel_pL(:,1:KaLim)./SfRel_vL(:,1:KaLim));               
Tsf_aR = abs(SfRel_pR(:,1:KaLim)./SfRel_vR(:,1:KaLim));

% Cm entry into unstable equilibrium
% Cm relative to Ba 
CUeRel_L = Cp_pL - Bp_pL;                       CUeRel_R = Cp_pR - Bp_pR;    % center of mass relative to stance foot

% CoM TtC - limited to first 50 frames to cut off end where high KaUE
% values arise (i.e. the jump is KaUE at the end of adaptation)
CmRel_pL = CUeRel_L - BaLhs_pL(:,1);            CmRel_pR = CUeRel_R - BaRhs_pR(:,1); % CoM relative to stance foot relative to Ba
CmRel_vL = [diff(CmRel_pL,[],2),diff(CmRel_pL(:,end-1:end),[],2)]./dt;  
CmRel_vR = [diff(CmRel_pR,[],2),diff(CmRel_pR(:,end-1:end),[],2)]./dt;

Tcm_aL = abs(CmRel_pL(:,1:KaLim)./CmRel_vL(:,1:KaLim));               
Tcm_aR = abs(CmRel_pR(:,1:KaLim)./CmRel_vR(:,1:KaLim));

% indexed position where CoM passes Bp
% obtaining zero crossings and differentiating to get 1's where cross occurs
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% accelCUrel_L = diff(CUeRel_L,[],2)*10; accelCUrel_L(isnan(accelCUrel_L)) = 0.001;
% accelCUrel_R = diff(CUeRel_R,[],2)*10; accelCUrel_R(isnan(accelCUrel_R)) = 0.001;
% 
% [~,cuLdiff] = max(diff(accelCUrel_L>0,[],2) == -1,[],2);    
% [~,cuRdiff] = max(diff(accelCUrel_R>0,[],2) == -1,[],2);    

% CfL = diff(CUeRel_L(:,1:cuLdiff)>0,[],2) == 1;      CfR = diff(CUeRel_R(:,1:cuRdiff)>0,[],2) == 1; 
CfL = diff(CUeRel_L>0,[],2) == 1;      CfR = diff(CUeRel_R>0,[],2) == 1; 

% Need to remove multiple 1's in each row
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%% Try removing the extra 1's on each row again
Cindex_L = [CfL,zeros(size(CfL,1),1)==1];     Cindex_R = [CfR,zeros(size(CfR,1),1)==1]; 

% Ka and Ka ratio at entry to UE
Ka_L = Tsf_aL./Tcm_aL;                          Ka_R = Tsf_aR./Tcm_aR;
KaUE_L = Ka_L(Cindex_L(:,1:KaLim));                        KaUE_R = Ka_R(Cindex_R(:,1:KaLim));
Ka_tL = Ka_L(:,1);                              Ka_tR = Ka_R(:,1);

% save data
dataKIN.Ka_tL = Ka_tL;                          dataKIN.Ka_tR = Ka_tR;
dataKIN.KaUE_L = KaUE_L;                        dataKIN.KaUE_R = KaUE_R;
%--------------------------------------------------------------------------

% plot saved parameters - Kp_t, Tcm, Tsf
if pltSav == 1
    subplot(2,2,3); hold on
    plot(1:length(Ka_tL),Ka_tL,'-k');              plot(1:length(Ka_tR),Ka_tR,'-r')
    xlabel('strides');                             ylabel({'anticipated boundary','ratio'})
    box on;                                        title('Ka_{t}')

    subplot(2,2,4); hold on
    plot(1:length(KaUE_L),KaUE_L,'-k');            plot(1:length(KaUE_R),KaUE_R,'-r');
    xlabel('strides');
    box on;                                        title('Ka_{UE}')
end

% plot check all parameters
if pltCHK == 1
    figure; hold on
    % ------------ plot Ka, Tcm, Tsf ------------------
    % left
    subplot(2,4,5); hold on
    yyaxis left;    p8 = plot(1:size(Ka_L,2),Ka_L,'-r');
    ylabel('Ka');   set(gca,'ycolor','k');  yy5 = get(gca,'ylim');
    ylim([0 min([yy5(2),4])]);              xx1 = get(gca,'xlim');
    if xx1(2) > 50; xlim([0 50]); end
    yyaxis right;   p9 = plot(1:size(Ka_L,2),Tcm_aL,'-k');
    p10 = plot(1:size(Ka_L,2),Tsf_aL,'-b');
    ylabel('TtC');  title('left limb');     ylim([0 1.5])
    legend([p8(1) p9(1) p10(1)],'Ka','Tcm','Tsf'); legend box off
    box on;         set(gca,'ycolor','k')
    xlabel('frames')
    
    % right
    subplot(2,4,6); hold on
    yyaxis left;    plot(1:size(Ka_R,2),Ka_R,'-r');
    ylabel('Ka');   set(gca,'ycolor','k'); yy6 = get(gca,'ylim');
    ylim([0 min([yy6(2),4])]);             xx2 = get(gca,'xlim');
    if xx2(2) > 50; xlim([0 50]); end
    yyaxis right;   plot(1:size(Ka_R,2),Tcm_aR,'-k');
    plot(1:size(Ka_R,2),Tsf_aR,'-b');
    ylabel('TtC');  title('right limb');    ylim([0 1.5])
    box on;         set(gca,'ycolor','k')
    xlabel('frames')
    
    % ------------ plot stride-by-stride KaUE,Kp_t,Tcm_t, Tsf_t ------------
    subplot(2,4,7); hold on
    yyaxis right
    p11 = plot(1:length(KaUE_L),KaUE_L,'-b');       plot(1:length(KaUE_R),KaUE_R,'--b');
    set(gca,'ycolor','b');                          ylabel('Ka @UE entry')
    yyaxis left
    p12 = plot(1:length(Ka_tL),Ka_tL,'-k');         plot(1:length(Ka_tR),Ka_tR,'--k')
    set(gca,'ycolor','k');      xlabel('strides');  ylabel('Ka @toe-off')
    legend([p11(1) p12(1)],'Ka @UE','Ka @Toe-off'); legend box off       
    box on;                                         title('Ka')
    
    subplot(2,4,8); hold on
    p12 = plot(1:length(Cm_crossL),Cm_crossL,'-k'); plot(1:length(Cm_crossR),Cm_crossR,'--k');
    p13 = plot(1:length(Sf_crossL),Sf_crossL,'-b'); plot(1:length(Sf_crossR),Sf_crossR,'--b');
    legend([p12(1) p13(1)],'Cm','Sf');    ylabel('% swing cross B_{p}'); xlabel('strides')
    box on;     ylim([0 100]); title('Cm and Sf @swing cross B_{p}'); legend box off
    
    
end

save(matfile,'dataKIN','-append')
end
