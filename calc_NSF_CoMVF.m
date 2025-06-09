function calc_NSF_CoMVF(matfile,V0,ii,plt) % [CoMvarY,CoMvarZ] = calc_CoMVar(file,fast)

clear dataCOM

load(matfile) 

[~,name,~] = fileparts(matfile);
eval(['Force = ' strrep(name,' ','_') '.Force;'])

fq = Force.Frequency;
dt = 1/fq;

LFx = dataForce.forceLeftX; % + dataForce.offsets.forceLeftX;
LFy = dataForce.forceLeftY; % + dataForce.offsets.forceLeftY;
LFz = dataForce.forceLeftZ; % + dataForce.offsets.forceLeftZ;

RFx = dataForce.forceRightX; % + dataForce.offsets.forceLeftX;
RFy = dataForce.forceRightY; % + dataForce.offsets.forceRightY;
RFz = dataForce.forceRightZ; % + dataForce.offsets.forceRightZ;

% Summed forces
Fx = LFx+RFx;
Fy = LFy+RFy;
Fz = LFz+RFz;

g = 9.8;

% Normalizing to weight
m = mean(Fz)/g;
% Calculating Acceleration
CoMX = (Fx)/m;
CoMY = (Fy)/m;
CoMZ = (Fz - m*g)/m;

% Getting contact events
LIC = Events.LICidx(1:end-1);
LTO = Events.LTOidx(1:end-1);
RIC = Events.RICidx(1:end-1);
RTO = Events.RTOidx(1:end-1);
% Trimming data to only within the full set of gait cycles to ensure
% periodic gait and summation to zero
[~,firstEvent] = min([LIC(1) LTO(1) RIC(1) RTO(1)]);
time = 1/fq:1/fq:size(dataForce.forceLeftX,1);
if firstEvent == 1
    cutFrames = LIC(1):LIC(end); Vt = time(LIC(1):LIC(end));
elseif firstEvent == 2
    cutFrames = LTO(1):LTO(end); Vt = time(LTO(1):LTO(end));
elseif firstEvent == 3
    cutFrames = RIC(1):RIC(end); Vt = time(RIC(1):RIC(end));
else
    cutFrames = RTO(1):RTO(end); Vt = time(RTO(1):RTO(end));
end

% lic = LIC(LIC<cutFrames(end)); lto = LTO(LTO<cutFrames(end));
% ric = RIC(RIC<cutFrames(end)); rto = RTO(RTO<cutFrames(end));

% Enforcing mean acceleration is zero
AoffsetZ = mean(CoMZ(cutFrames));
AoffsetY = mean(CoMY(cutFrames));
AoffsetX = mean(CoMX(cutFrames));

% Enforcing mean vertical velocity is zero
vZeroZ = mean(cumtrapz(CoMZ(cutFrames) - AoffsetZ)*dt);
vZeroY = mean(cumtrapz(CoMY(cutFrames) - AoffsetY)*dt + V0);
vZeroX = mean(cumtrapz(CoMX(cutFrames) - AoffsetX)*dt);

% Calculating COM velocity
dataCOM.Vz = cumtrapz(CoMZ(cutFrames) - AoffsetZ)*dt - vZeroZ;
dataCOM.Vy = cumtrapz(CoMY(cutFrames) - AoffsetY)*dt + V0 - vZeroY;
dataCOM.Vx = cumtrapz(CoMX(cutFrames) - AoffsetX)*dt - vZeroX;

lic = LIC(LIC<length(dataCOM.Vz))-cutFrames(1)+1; 
lto = LTO(LTO<length(dataCOM.Vz))-cutFrames(1)+1;
ric = RIC(RIC<length(dataCOM.Vz))-cutFrames(1)+1; 
rto = RTO(RTO<length(dataCOM.Vz))-cutFrames(1)+1;

dataCOM.Eventoffset = cutFrames(1);

dataCOM.CoMZ = CoMZ(cutFrames); 
dataCOM.CoMY = CoMY(cutFrames);
dataCOM.CoMX = CoMX(cutFrames);

dataCOM.LFz = LFz(cutFrames)/1000;
dataCOM.LFy = LFy(cutFrames)/1000;
dataCOM.LFx = LFx(cutFrames)/1000;
dataCOM.RFz = RFz(cutFrames)/1000;
dataCOM.RFy = RFy(cutFrames)/1000;
dataCOM.RFx = RFx(cutFrames)/1000;


if plt == 1
figure(ii); clf; 
subplot(211); hold on
plot(dataCOM.Vz,'b')
plot(lic,dataCOM.Vz(lic),'ok')
plot(lto,dataCOM.Vz(lto),'xk')
plot(ric,dataCOM.Vz(ric),'or')
plot(rto,dataCOM.Vz(rto),'xr')

plot(LFz(cutFrames)/1000,'k')
plot(RFz(cutFrames)/1000,'r')

subplot(212); hold on
plot(dataCOM.CoMZ,'b')
plot(lic,dataCOM.CoMZ(lic),'ok')
plot(lto,dataCOM.CoMZ(lto),'xk')
plot(ric,dataCOM.CoMZ(ric),'or')
plot(rto,dataCOM.CoMZ(rto),'xr')

% figure(2); clf; hold on
% plot(dataCOM.Vy,dataCOM.Vz,'b')
% plot(dataCOM.Vy(LIC),dataCOM.Vz(LIC),'ok')
% plot(dataCOM.Vy(LTO),dataCOM.Vz(LTO),'xk')
% plot(dataCOM.Vy(RIC),dataCOM.Vz(RIC),'or')
% plot(dataCOM.Vy(RTO),dataCOM.Vz(RTO),'xr')
end

save(matfile,'dataCOM','-append')
end
        
