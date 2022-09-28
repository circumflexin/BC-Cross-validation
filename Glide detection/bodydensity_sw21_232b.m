% % % % %%%Body density estimation. Japan workshop May 2016
% % % %Lucia Martina Martin Lopez lmml2@st-andrews.ac.uk
% % % % Edited by Eilidh Siegal, July 2017 (es250@st-andrews.ac.uk)
% % % % All edits by ES are commented, with comments preceded by "ES - "
% % % 
% % % % Edited by Alec Burslem, 2019-2024 (acb35@st-andrews.ac.uk)
% % % AB - Hardcoded accelerometer/magnetometer axes for cetaceans - to use for seals, change the dominant axis in magnet_rot_sa() and Ahf_Anlf() calls.
% % % AB - Added calculation of glide/fluke cutoff from terminal glides
% % % AB - Added diagnostic plots
% % % AB - Added buffer at start of detected glides
% % % AB - New splitGL and SWdensityFromCTD functions 
% % % AB - Use Magnetometer method to find fluke rate
% % % AB - Exclude foraging periods from ascent and descent phase to make glide ratio more accurate
% % % Ab - Reordered so manual dive phase classification is done first 
% % % AB - Version controlled repo on github: @circumflexin
% % % AB - new magnet_rot_AB to waste fewer gliding periods

clc                      
clear           
%% Config
%LAT = % Latitude where CTD measurements took place
Glide_Sum_Dir = 'D:\Data\BC_crossval\Glide_detection\Summaries\';
Glide_Ratio_Dir =  'D:\Data\BC_crossval\Glide_detection\Ratios';
alpha=25; %degrees, cuts off magnetometry when the sensor is too closely aligned with the earth's magnetic field

%% LOAD DATA
% Load prh
tag='sw21_232b';                 % Insert deployment here
settagpath('prh','D:\Raw\Azores_2021\MDTAGs\Data\D3\metadata\prh')
loadprh(tag);
%Load CTD
%CTD = readtable(''); % 
% inspect OTAB
settagpath('cal','D:\Raw\Azores_2021\MDTAGs\Data\D3\metadata\cals')
[CAL,DEPLOY,ufname] = d3loadcal(tag)
DEPLOY.OTAB

%% DEFINE DIVES 
% and make a summary table describing the characteristics of each dive.

% ES - used 2m as the threshold of the dive (surface) and 10m for dive definition (i.e. deep dive if below 10m)
% as in the IgorPro FindDives.ipf provided by Tomoko Narazaki (edited last 19 May 2016)

k=1:length(p);
mindivedef = 100;                       % AB - to exclude resting dives; where bouyancy is actively maintained. 
surface=7.5;                            % AB - increased due to longer body
T=finddives(p,fs,mindivedef,surface); 
D=[];                                   % [start_time(s) end_time(s) dive_duration(s) max_depth(s) max_depth(m) ID(n)]
D(:,1)=T(:,1);                          % start time of each dive in seconds since the tag on time
D(:,2)=T(:,2);                          % end time of each dive in seconds since the tag on time
D(:,3)=T(:,2)-T(:,1);                   % dive duration in seconds
D(:,4)=[T(2:end,1)-T(1:end-1,2);NaN];   % post-dive surface duration in seconds
D(:,5)=T(:,4);                          % time of the deepest point of each dive
D(:,6)=T(:,3);                          % maximum dive depth of each dive
D(:,7)=1*(1:size(T,1))';                % dive ID number, starting from 1
clear T

%4_DEFINE PRECISE DESCENT AND ASCENT PHASES 
Bottom=[];
Phase(1:length(p))=NaN;
DES=[];
ASC=[];

isdive = false(size(p));
for i = 1:size(D,1); isdive(round(D(i,1)*fs):round(D(i,2)*fs)) = true; end
pd = p; pnd = p; pd(~isdive) = nan; pnd(isdive) = nan;
I = 1:length(p);
figure(4); clf; 
sp1 = subplot(5,1,1:4)
plot(I,pd,'g'); hold on; plot(I,pnd,'r'); set(gca,'ydir','rev','ylim',[min(p) max(p)]);
tit = title('Click within the last dive you want to use');
legend('Dive','Not Dive','location','southeast');
sp2=subplot(5,1,5)
plot(I,pitch,'g');
legend('pitch')
x = ginput(1);
nnn = find(D(:,1)<x(1)/fs,1,'last');

for dive=1:nnn                                 % ES - Line added so as to only include until dive where clicked previously on dive profile
    kk=round(fs*D(dive,1)):round(fs*D(dive,2)); % it is selecting the whole dive
    %kkI = false(size(p)); kkI(kk) = true;
      %enddes=((find((smoothpitch(kk)*180/pi)>0,1,'first')+D(dive,1)*fs));% search for the first point at which pitch is positive
      %startasc=((find((smoothpitch(kk)*180/pi)<0,1,'last')+D(dive,1)*fs));%search for the last point at which the pitch is negative
    %    if you want to do it manually as some times there is a small ascent
    %     phase during the descent and a small descent phase during the ascent.
    figure(5)
    ax(1)=subplot(211); plott(p(kk),fs)	% plott plots sensor data against a time axis
    title(join("Dive: "+dive))
    ax(2)=subplot(212); plott(pitch(kk)*180/pi,fs,0)
    linkaxes(ax, 'x'); % links x axes of the subplots for zoom/pan
    [x,y]=ginput(2);	% click on where the pitch angle first goes to zero in the descent and last goes to zero in the ascent
    enddes=round(x(1))*fs+D(dive,1)*fs;         % ES - changed output names and multiples by fs to get correct time frame
    startasc=round(x(2))*fs+D(dive,1)*fs;       % ES - changed output names and multiples by fs to get correct time frame
    Phase(kk(kk<enddes))=-1; %Descent
    Phase(kk(kk<startasc & kk>enddes)) = 0 ; % Bottom
    Phase(kk(kk>startasc)) = 1 ; % Ascent
    Bottom(dive,1)=(enddes)/fs; %Time in seconds at the start of bottom phase (end of descent)
    Bottom(dive,2)= p(enddes);% Depth in m at the start of the bottom phase (end of descent phase)
    Bottom(dive,3)=(startasc)/fs;%Time in seconds at the end of bottom phase (start of ascent)
    Bottom(dive,4)= p(startasc);% Depth in m at the end of the bottom phase (start of ascent phase)
    des=round(fs*D(dive,1)):round(enddes);%selects the whole descent phase
    asc=round(startasc):round(fs*D(dive,2));%selects the whole ascent phase
    DES=[DES,des];                                               % Concats all descents
    ASC=[ASC,asc];                                               % Concats all descents
end
pasc = p; pdes = p;
pasc(Phase<1 | isnan(Phase)) = nan;
pdes(Phase>-1 | isnan(Phase)) = nan;
plot(sp1,pasc,'k'); plot(sp1,pdes,'b');
legend('Dive','Not Dive','Ascent','Descent');
delete(tit)
linkaxes([sp1 sp2], 'x');

dives_w_exclusions = [3]

% AB - exclude any periods? This lets you exclude foraging periods on the
% way up/down from glude:stroke ratio analyses
 for n = dives_w_exclusions % list dive IDs you want to exclude sections from, to exclude multiple sections, enter the dive twice
     kk=round(fs*D(n,1)):round(fs*D(n,2));
     figure(6)
     ax(1)=subplot(211); plott(p(kk),fs)	% plott plots sensor data against a time axis
     ax(2)=subplot(212); plott(pitch(kk)*180/pi,fs,0)
     linkaxes(ax, 'x'); % links x axes of the subplots for zoom/pan
     [x,y]=ginput(2);	% click on start and end of period to exclude
     start_exc = round(x(1))*fs+D(n,1)*fs;
     end_exc = round(x(2))*fs+D(n,1)*fs;
  
   Phase(kk(kk<end_exc & kk>start_exc)) = 0;
 end
isdive = false(size(p));
for i = 1:size(D,1); isdive(round(D(i,1)*fs):round(D(i,2)*fs)) = true; end
pd = p; pnd = p; pd(~isdive) = nan; pnd(isdive) = nan;
I = 1:length(p);
figure(4); clf; 
sp1 = subplot(5,1,1:3);
plot(I,pd,'g'); hold on; plot(I,pnd,'r'); set(gca,'ydir','rev','ylim',[min(p) max(p)]);
tit = title('Click within the last dive you want to use');
legend('Dive','Not Dive','location','southeast');
hold on
pasc = p; pdes = p;
pasc(Phase<1 | isnan(Phase)) = nan;
pdes(Phase>-1 | isnan(Phase)) = nan;
plot(pasc,'k'); plot(pdes,'b');
legend('Dive','Not Dive','Ascent','Descent','location','southeast');
delete(tit)

%% save asc/des cues here to avoid having to repeat manual classification
%save(fullfile("D:\Data\BC_crossval\Glide_detection\", tag))
tag='sw21_232b';   
load(fullfile("D:\Data\BC_crossval\Glide_detection\", tag))

%% 3.2_SEPARATE LOW AND HIGH ACCELERATION SIGNALS 

%3.2.1_select periods of analysis.
% The power spectra is calculated of the longitudinal and 
% dorso-ventral accelerometer signals during descents and ascents to 
% determine the dominant stroke frequency for each animal in each phase with
% a fft of 512 and a sampling rate of fs. 
% Output: S is the amount of power in each particular frequency (f)

% During whole deployment 
[S,f]=speclev(Aw(k,:),512,fs);plot(f,S);
% During all descents and ascents phases where mainly steady swimming
% occurs. When calculated for the whole dive it may be difficult to 
% differenciate the peak at which stroking rate occurs as there is other 
% kind of movements than only steady swimming

% Using magnetometer method
[~,~,Sa] = magnet_rot_sa_AB(Aw,Mw,fs,0.05,alpha,1,1:length(p),[],[]); %fixed cutoff
Sa(isnan(Sa))=0; % Convert Nans into zeroes or it breaks speclev
[S,f]=speclev(Sa(k,:),512,fs);
plot(f,S)

ACCP=Aw; % Sa or Aw

%3.2.2_Determine FR fluking rate and cut-off frequency
% calculate the power spectrum of the accelerometer data at the whale frame
% plot PSD for the whole deployment, all descents and all ascent phases.
fluke_rate = nan(6,1); FRmag = fluke_rate; 
FRl = nan(6,1); FRlmag = FRl;
FRn = 1;
cs = 'bgr';
va = {'top';'';'bottom'};

% Subplot 1: whole deployment
figure (1); clf; 
[S,f]=speclev(ACCP(k,:),512,fs);

ax(1)=subplot(311);
for nn = [1 3]
    [peakloc, peakmag] = peakfinder(S(:,nn),0.1); % make 0.1 smaller if it is not finding the peak you expect
    peakloc(1) = []; peakmag(1) = [];
    smoothS = runmean(S(:,nn),5);
    peakmag = peakmag - smoothS(peakloc);
    [~,peak] = max(peakmag);
    peak = peakloc(peak);
    plot(f(peak),S(peak,nn),'go','markersize',10,'linewidth',2); 
    [minf,bf] = min(S(1:peak,nn));
    FRl(FRn) = f(bf);
    hold on
    plot([f(bf) f(bf)],[min(S(:,nn)) minf],'k--','linewidth',2);
    text(f(bf),min(min(S(:,[1 3]))),['f = ' num2str(round(f(bf)*100)/100)],'horizontalalignment','right','color',cs(nn),'verticalalignment',va{nn});
    fluke_rate(FRn) = f(peak);
    FRmag(FRn) = S(peak,nn);
    FRn = FRn+1;
end
[~,b] = max(FRmag(1:2));
if b == 2; v1 = 'top'; v2 = 'bottom'; else v1 = 'bottom'; v2 = 'top'; end
text(fluke_rate(1),FRmag(1)+2*sign(length(v1)-length(v2)),num2str(round(fluke_rate(1)*100)/100),'verticalalignment',v1,'horizontalalignment','center');
text(fluke_rate(2),FRmag(2)-2*sign(length(v1)-length(v2)),num2str(round(fluke_rate(2)*100)/100),'verticalalignment',v2,'horizontalalignment','center');
b = plot(f,S(:,1),'b'); grid;
r = plot(f,S(:,3),'r');
grid on;
ax1=gca;
set(get(ax1,'Xlabel'),'String', [{'$\bf\ Frequency \hspace{1mm} (Hz) $'}],'interpreter','latex','FontSize',8,'FontName','Arial')
ys = get(ax1,'ylim');
hp1 = text(0.02,diff(ys)*.92+min(ys),'Whole deployment','FontSize',12,'FontWeight','bold','horizontalalignment','left') ;% to write the name of the pannel
legend([b r],'HPF acc x axis (surge)','HPF acc z axis (heave)')

% Subplot 2: descent
ax(2)=subplot(312);
[S,f]=speclev(ACCP(DES,:),512,fs);
for nn = [1 3]
    [peakloc, peakmag] = peakfinder(S(:,nn),0.001); % make 0.1 smaller if it is not finding the peak you expect
    peakloc(1) = []; peakmag(1) = [];
    smoothS = runmean(S(:,nn),10);
    peakmag = peakmag - smoothS(peakloc);
    [~,peak] = max(peakmag);
    peak = peakloc(peak);
    %plot(f(peak),S(peak,nn),'go','markersize',10,'linewidth',2);
    %[minf,bf] = min(S(1:peak,nn));
    %FRl(FRn) = f(bf);
    hold on
    plot([f(bf) f(bf)],[min(S(:,nn)) minf],'k--','linewidth',2);
    text(f(bf),min(min(S(:,[1 3]))),['f = ' num2str(round(f(bf)*100)/100)],'horizontalalignment','right','color',cs(nn),'verticalalignment',va{nn});
    fluke_rate(FRn) = f(peak);
    FRmag(FRn) = S(peak,nn);
    FRn = FRn+1;
end
[~,b] = max(FRmag(1:2));
if b == 2; v1 = 'top'; v2 = 'bottom'; else v1 = 'bottom'; v2 = 'top'; end
text(fluke_rate(1),FRmag(3)+2*sign(length(v1)-length(v2)),num2str(round(fluke_rate(3)*100)/100),'verticalalignment',v1,'horizontalalignment','center');
text(fluke_rate(2),FRmag(4)-2*sign(length(v1)-length(v2)),num2str(round(fluke_rate(4)*100)/100),'verticalalignment',v2,'horizontalalignment','center');
b = plot(f,S(:,1),'b'); grid;
r = plot(f,S(:,3),'r');
grid on;
ax1=gca;
set(get(ax1,'Xlabel'),'String', [{'$\bf\ Frequency \hspace{1mm} (Hz) $'}],'interpreter','latex','FontSize',8,'FontName','Arial')
ys = get(ax1,'ylim');
hp1 = text(0.02,diff(ys)*.92+min(ys),'Descents','FontSize',12,'FontWeight','bold','horizontalalignment','left') ;% to write the name of the pannel
legend([b r],'HPF acc x axis (surge)','HPF acc z axis (heave)')

% Subplot 3: ascent
ax(3)=subplot(313);
[S,f]=speclev(ACCP(ASC,:),512,fs);
for nn = [1 3]
[peakloc, peakmag] = peakfinder(S(:,nn),0.001);
peakloc(1) = []; peakmag(1) = [];
smoothS = runmean(S(:,nn),10);
peakmag = peakmag - smoothS(peakloc);
[~,peak] = max(peakmag);
peak = peakloc(peak);
plot(f(peak),S(peak,nn),'go','markersize',10,'linewidth',2)
[minf,bf] = min(S(1:peak,nn));
FRl(FRn) = f(bf);
hold on
plot([f(bf) f(bf)],[min(S(:,nn)) minf],'k--','linewidth',2);
text(f(bf),min(min(S(:,[1 3]))),['f = ' num2str(round(f(bf)*100)/100)],'horizontalalignment','right','color',cs(nn),'verticalalignment',va{nn});
fluke_rate(FRn) = f(peak);
FRmag(FRn) = S(peak,nn);
FRn = FRn+1;
end
[~,b] = max(FRmag(3:4));
if b == 2; v1 = 'top'; v2 = 'bottom'; else v1 = 'bottom'; v2 = 'top'; end
text(fluke_rate(3),FRmag(5)+2*sign(length(v1)-length(v2)),num2str(round(fluke_rate(5)*100)/100),'verticalalignment',v1,'horizontalalignment','center');   
text(fluke_rate(4),FRmag(6)-2*sign(length(v1)-length(v2)),num2str(round(fluke_rate(6)*100)/100),'verticalalignment',v2,'horizontalalignment','center'); 
b = plot(f,S(:,1),'b'); grid;
r = plot(f,S(:,3),'r');
grid on;
ax2=gca;
set(get(ax2,'Xlabel'),'String', [{'$\bf\ Frequency \hspace{1mm} (Hz) $'}],'interpreter','latex','FontSize',8,'FontName','Arial')
ys = get(ax2,'ylim');
hp1 = text(0.02,diff(ys)*.92+min(ys),'Ascents','FontSize',12,'FontWeight','bold','horizontalalignment','left') ;% to write the name of the pannel
legend([b r],'HPF acc x axis (surge)','HPF acc z axis (heave)')
linkaxes(ax, 'x'); % links x axes of the subplots for zoom/pan




% f = number that multiplied by the FR gives the cut-off frequency fl, of 
% the low pass filter. f is a fraction of FR. You can set default value to
% 0.4 or set f as fl (frequency at the negative peak in the
% power spectral density plot)/FR.
fluke_rate = mean(fluke_rate([4,5])); % Change indexing to calculate from a subset of the data (i.e only ascents)
cutoff = mean(0.09) % Change indexing to calculate from a subset of the data (i.e only ascents)
%or enter manually if the peak detection isn't working
%fluke_rate = mean([0.1953,0.1953]);  % Change indexing to calculate from a subset of the data (i.e only ascents)
f = cutoff/fluke_rate;

%transition: 0.078 - 0.09

%3.2.3_Separate low and high pass filtered acceleration signal using the
%parameters defined earlier and the function Ahf_Alnf.
[Anlf,Ahf,GL,KK] = Ahf_Anlf(Aw,fs,fluke_rate,f,3,k,[],[]) ;

%3.2.4_ calculate the smooth pitch from the low pass filter acceleration signal
% to avoid incorporating signals above the stroking frequency
[smoothpitch,smoothroll]=a2pr(Anlf(k,:));
%check the difference between pitch and smoothpitch
figure (2); clf;
ax1 = subplot(2,1,1);
plott(p,fs)
ax2 = subplot(2,1,2);
plot((1:length(p))/fs,pitch*180/pi)
hold on
plot((1:length(p))/fs,smoothpitch*180/pi,'r')
legend('pitch','smoothpitch')
linkaxes([ax1 ax2], 'x');

%% 5_ESTIMATE SWIM SPEED

thdeg = 30; %degree threshold above which speed can be estimated % increase for sperm whales?
[SwimSp] = inst_speed_AB(p,smoothpitch,fs,cutoff,k,thdeg); % AB - changed pitch to sin(pitch), takes cutoff directly
%help inst_speed
figure(4);
sp2= subplot(5,1,5);
plot(k,SwimSp,'b'); ylim([0 max(SwimSp)]);
legend('speed')
linkaxes([sp1 sp2], 'x'); % links x axes of the subplots for zoom/pan

thdeg = 30; %degree threshold above which speed can be estimated % increase for sperm whales?
[SwimSp] = inst_speed_AB(p,smoothpitch,fs,cutoff,k,thdeg); % AB - changed pitch to sin(pitch), takes cutoff directly
%help inst_speed

%inspect in 
figure(4); clf; 
sp1 = subplot(5,1,1:2);
plot(I,pd,'g'); hold on; plot(I,pnd,'r'); set(gca,'ydir','rev','ylim',[min(p) max(p)]);
tit = title('Click within the last dive you want to use');
legend('Dive','Not Dive','location','southeast');
hold on
pasc = p; pdes = p;
pasc(Phase<1 | isnan(Phase)) = nan;
pdes(Phase>-1 | isnan(Phase)) = nan;
plot(pasc,'k'); plot(pdes,'b');
legend('Dive','Not Dive','Ascent','Descent','location','southeast');
delete(tit)
sp2 = subplot(5,1,3);
plot((1:length(p)),pitch*180/pi)
hold on
plot((1:length(p)),smoothpitch*180/pi,'r')
legend('pitch','smoothpitch')
hold on 
sp3= subplot(5,1,4:5);
plot(k,SwimSp,'b'); ylim([0 max(SwimSp)]);
legend('speed')
linkaxes([sp1 sp2 sp3], 'x'); % links x axes of the subplots for zoom/pan


%% 6_ESTIMATE SEAWATER DESNSITY AROUND THE TAGGED WHALES

% % ES - edited script here to include pressure (i.e. depth) in the SWdensity estimation.
% % Involved changing sw_dens0(sali,temp) to sw_dens(sali,temp,P). 
% % 1st step required calculating pressure from depth (because pressure, not depth, is input into sw_dens fucntion)
% % Downloaded gsw_matlab_v3_05_7 tools from: http://www.teos-10.org/software.htm
% % Used gsw_p_from_z function (edited to gsw_p_from_z_ES) to get pressure from depth
% % Downloaded Mixing Oceanographic toolbbox from: https://uk.mathworks.com/matlabcentral/fileexchange/47595-mixing--mx--oceanographic-toolbox-for-em-apex-float-data/content/mixing_library/private1/seawater/sw_dens.m
% % Used sw_dens tool within SWdesnityfromCTD_ES to calculate SW density with pressure (i.e. depth) considered
% % lat = 70.82611; latitude of CTD taken from Miller et al 2016 and  converted to decimal degrees needed for gsw_matlab_v3_05_7
% 
% % Home > Import > CTD data as matrix
% [val index] = max(CTD.Press)
% CTDu = CTD(index:height(CTD),:) %upcast
% CTDd = CTD(1:index,:) % downcast
% CTDf = CTDu
% 
% PRE = CTDf{:,'Press'};
% DPT = pres2depth(CTDf{:,'Press'}, LAT); % AB -  Assumes decibars. Convert CTD data if not, uses Paul Wensveen's function to calculate depth.
% 
% %Seawater density estimation
% TMP = CTDf{:,'Temp'}; % AB
% SL = CTDf{:,'Sal_'}; %  AB
% %SWdensity1 = CTD{:,'Density'};
% 
% 
% 
% [SWdensity,depCTD,temp]=SWdensityFromCTD_AB(DPT, TMP, SL,D, LAT, PRE);        % ES - edited function to make temp an output and to include depth (pressure) in seawater density calculation AB - Edited to accept sverdrup CTD data and use the raw values for pressure.
% Dsw=EstimateDsw(SWdensity, depCTD, p);
% tempw = Estimate_temp_ES(temp, depCTD, p);                          % ES - added line and created function to calculate temperature around the whale


%% 7_EXTRACT STROKES AND GLIDES
% Can be done using the body rotations (pry) estimated using the
% magnetometer method, or it can be done using the dorso-ventral axis of the
% high-pass filtered acceleration signal
% Using whichever method, tmax and J need to be determined.
% tmax is the maximum duration allowable for a fluke stroke in seconds, it can be set as 1/FR
% J is the magnitude threshold for detecting a fluke stroke in m /s^2


%7.1  Using the heave high pass filtered acceleration signal,(use n=3)
% units of Ja are in m/s2; set Ja and tmax [] until determined.
[~,Ahf,~,~] = Ahf_Anlf(Aw,fs,fluke_rate,f,3,k,[],[]);

figure (5); clf;
sp1 = subplot(1,4,1:3);
ax2 = plotyy(1:length(p),p,1:length(pitch),Ahf(:,3));
set(ax2(1),'ydir','rev','ylim',[0 max(p)]);
set(ax2(2),'nextplot','add');
plot(ax2(2),1:length(pitch),Ahf(:,1),'m');
maxy = max(max(Ahf(round(D(1,1)*fs):round(D(nn,2)*fs),[1 3])));
set(ax2(2),'ylim',[-2*maxy 2*maxy],'ytick',round(-2*maxy*10)/10:0.1:2*maxy)
legend('Depth','HPF acc z axis','HPF acc x axis');
title('Zoom in to find appropriate thresholds for fluking, then enter it for J');
linkaxes([ax1 ax2], 'x');

flukeAcc1A = Ahf(Phase == 1,1); % hpf-x acc ascent
flukeAcc1D = Ahf(Phase == -1,1); % hpf-x acc descent
flukeAcc3A = Ahf(Phase == 1,3); % hpf-z acc ascent
flukeAcc3D = Ahf(Phase == -1,3); % hpf-z acc descent

sp2 = subplot(2,4,4);
TOTAL = abs([flukeAcc1A; flukeAcc1D]);  % to look at heave, change to 3
Y=buffer(TOTAL,2*round(1/fluke_rate*fs));
hist(max(Y),100); set(sp2,'xlim',[0 max(max(Y))]);
title('hpf-x acc');

sp3 = subplot(2,4,8);
TOTAL = abs([flukeAcc3A; flukeAcc3D]);
Y=buffer(TOTAL,2*round(1/fluke_rate*fs));
hist(max(Y),100,'FaceColor','g'); set(sp3,'xlim',[0 max(max(Y))]);
title('hpf-z acc');


% Choose a value for J based on the histogram for:
%   hpf-x, then when detecting glides in the next step use Ahf_Anlf
%   function with n=1
%   hpf-z then when detecting glides in the next step use Ahf_Anlf
%   funct ion with n=3

Ja= 0.09;    % in m/s2
tmax=1/fluke_rate;  % in seconds
[Anlf,Ahf,GLa,KKa] =  Ahf_Anlf(Aw,fs,fluke_rate,f,3,k,Ja,tmax);
GLa(:,1) = GLa(:,1) + (1/fluke_rate)/4;



%% 7.2 Using the body rotations (pry) estimated using the magnetometer method (use n=1)
figure (6); clf;
sp1 = subplot(1,4,1:3);
[MagAcc,pry,Sa] = magnet_rot_sa(Aw,Mw,fs,fluke_rate,f,alpha,1,k,[],[]);
ax1 = plotyy(1:length(p),p,1:length(pry(:,1)),pry(:,1));
set(ax1(1),'ydir','rev','ylim',[0 max(p)]);
maxy = max(pry(round(D(1,1)*fs):round(D(nn,2)*fs)));
set(ax1(2),'ylim',[-2*maxy 2*maxy],'ytick',round(-2*maxy*10)/10:0.1:2*maxy)
set(ax1(2),'nextplot','add');
legend('Depth','rotations in y axis');
title('Zoom in to find appropriate thresholds for fluking, then enter it for J');
sp2 = subplot(1,4,4);
flukePA = pry(Phase == 1,1); % body rotations ascent in radians
flukePD =pry(Phase == -1,1);
TOTAL = abs([flukePA; flukePD]);
Y=buffer(TOTAL,2*round(1/fluke_rate*fs));
hist(max(Y),100); set(sp2,'xlim',[0 max(max(Y))]);
% Choose a value for J based on the histogram for:
%   pry(:,1), then when detecting glides in the next step use magnet_rot_sa
%   function with n=1

%%
%Dev terminal glide detection: use this to calculate the body rotation
%during terminal glides, as a potential minimum value for J.

term_glide = [0,0]
term_glide_J_max = []

figure
for dive=1:nnn
AB_asc = round(fs*Bottom(dive,3)):round(fs*D(dive,2));%selects the whole ascent phase
sp1=subplot(211); plott(p(AB_asc),fs);	% plott plots sensor data against a time axis
pry_asc = pry(AB_asc);
t=(0:length(pry(AB_asc))-1)/fs;
sp2=subplot(212); plot(t,pry(AB_asc));
linkaxes([sp1,sp2], 'x'); % links x axes;
term_glide=ginput(2);
term_glide(:,1) = term_glide(:,1)*fs;
term_glide_J = pry_asc(term_glide(1,1):term_glide(2,1));
term_glide_J = term_glide_J.^2;
peaks = peakfinder(term_glide_J);
term_glide_J_max(dive) = sqrt(max(term_glide_J));
hold on;
plot((peaks+term_glide(1,1)-1)/fs,pry_asc(round(peaks+term_glide(1,1))),'y*');
waitforbuttonpress
hold off
end

hold off
figure;
plot(term_glide_J_max)
median(term_glide_J_max)
tmax=1/fluke_rate;%in seconds


% if the tag slips 
% %Jm = 0.0292;% in radians
% 
% Jm1 = 0.005;% in radians
% Aw1 = Aw(1:round(DEPLOY.OTAB(2,1)*fs),:);
% Mw1 = Mw(1:round(DEPLOY.OTAB(2,1)*fs),:);
% k1 = k(1:round(DEPLOY.OTAB(2,1)*fs));
% 
% Jm2 = 0.005;% in radians
% Aw2 = Aw(round(DEPLOY.OTAB(2,2)*fs):length(Aw),:);
% Mw2 = Mw(round(DEPLOY.OTAB(2,2)*fs):length(Aw),:);
% k2 = k(round(DEPLOY.OTAB(2,2)*fs):length(Aw));
% 
% [~,pry,~,~,~] = magnet_rot_sa_AB(Aw,Mw,fs,cutoff,alpha,1,k,[],[]);
% [~,~,~,GLm1,KKm1] = magnet_rot_sa_AB(Aw1,Mw1,fs,cutoff,alpha,1,k1,Jm1,tmax);
% [~,~,~,GLm2,KKm2] = magnet_rot_sa_AB(Aw,Mw,fs,cutoff,alpha,1,k2,Jm2,tmax);
% KKm2 = KKm2+DEPLOY.OTAB(2,2)
% KKm = [KKm1;KKm2]
% GLm2 = GLm2+DEPLOY.OTAB(2,2)
% GLm = [GLm1;GLm2]
% 


Jm = 0.0059;% in radians
tmax=1/fluke_rate;%in seconds
[~,pry,~,GLm,KKm] = magnet_rot_sa_AB(Aw,Mw,fs,cutoff,alpha,1,k,Jm,tmax,true);
GLm(:,1) = GLm(:,1) + (1/fluke_rate)/4; % AB - add half the fluking period to the glide start time, does not produce negatives because glides already have to be > fluke period.
%any((sqrt((GLm(:,2)-GL(:,1)).^2) == GLm(:,2)-GLm(:,1))==false) % check for negatives


%KKm/a = matrix of cues to zero crossings in seconds (1st column) and
%zero-crossing directions (2nd column). +1 means a positive-going
%zero-crossing. Times are in seconds.
%this is already ensuring that all glides are longer than tmax/2

%Check glides duration and positive and negative zero crossings (KK) based
% on selected J and tmax
% magnetometer method
figure (7); clf;
sp1 = subplot(5,1,1:4);
plot(k,pd,'g'); hold on; plot(k,pnd,'r','linewidth',2); set(gca,'ydir','rev','ylim',[min(p) max(p)]);
legend('Dive','Not Dive','location','southeast');
plot(sp1,pasc,'k','linewidth',2); plot(sp1,pdes,'b','linewidth',2);


sp2= subplot(5,1,5);
ax2= plotyy(k,SwimSp,k,pry(:,1));
set(ax2(2),'nextplot','add')
plot(ax2(2),KKm(:,1)*fs,pry(round(KKm(:,1)*fs),1),'r*')
set(ax2(1),'ylim',[0 max(SwimSp)]);
set(ax2(2),'ylim',[min(pry(Phase==-1,1)), max(pry(Phase==1,1))]);
legend('speed','body rotations')
linkaxes([sp1 sp2 ax2], 'x'); % links x axes

GLdurm= GLm(:,2)-GLm (:,1);
GLTm=[GLm(:,1),GLdurm];
GLdurm= GLm(:,2)-GLm (:,1);
GLTm=[GLm(:,1),GLdurm];
t=(0:length(pry(:,1))-1)/fs;
glkm=eventon(GLTm,t);
pgl=p;
pgl(glkm==0)=NaN;
h2=plot(sp1,t*fs,pgl,'-y','Linewidth',3);% glides detected with the body rotations (pry)
hold on
legend(sp1,'Bottom','Not Dive','Ascent','Descent','Glide','location','southeast');

% Accleromter method
figure (8); clf;
sp1 = subplot(5,1,1:4);
plot(k,pd,'g'); hold on; plot(k,pnd,'r','linewidth',2); set(gca,'ydir','rev','ylim',[min(p) max(p)]);
legend('Dive','Not Dive','location','southeast');
plot(sp1,pasc,'k','linewidth',2); plot(sp1,pdes,'b','linewidth',2);

sp2= subplot(5,1,5);
ax2= plotyy(k,SwimSp,k,Ahf(:,3));
set(ax2(2),'nextplot','add')
%plot(ax2(2),KKa(:,1)*fs,0,'r*')
yline(ax2(2),Ja,'r','linewidth',1);
yline(ax2(2),-Ja,'r','linewidth',1);
set(ax2(1),'ylim',[0 max(SwimSp)]);
set(ax2(2),'ylim',[min(Ahf(Phase==-1,1)), max(Ahf(Phase==1,1))]);
legend('speed','Highpass Acceleration')
linkaxes([sp1 sp2 ax2], 'x'); % links x axes

GLdura= GLa(:,2)-GLa (:,1);
GLTa=[GLa(:,1),GLdura];
t=(0:length(pry(:,1))-1)/fs;
glka=eventon(GLTa,t);
pgl=p;
pgl(glka==0)=NaN;
h3=plot(sp1,t*fs,pgl,'-y','Linewidth',3);%
legend(sp1,'Bottom','Not Dive','Ascent','Descent','Glide','location','southeast');

% Comparison 
agree = glka+glkm;
agree_yes = p; agree_no = p;onlya = p; onlym = p;
agree_yes(agree~=2)=NaN;
onlya(agree~=1|glka==0)=NaN;
onlym(agree~=1|glkm==0)=NaN;

figure (9); clf;
sp1 = subplot(6,1,1:4);
plot(k,p,'k','Linewidth',3) ; hold on; set(gca,'ydir','rev','ylim',[min(p) max(p)]);
hold on
plot(t*fs,agree_yes,'-y','Linewidth',3);
plot(t*fs,onlym,'-m','Linewidth',3);
plot(t*fs,onlya,'-r','Linewidth',3);
legend('Track','Agree','Only Rot' ,'Only Acc');

sp2= subplot(6,1,5);
ax2= plotyy(k,SwimSp,k,pry(:,1));
set(ax2(2),'nextplot','add')
plot(ax2(2),KKm(:,1)*fs,pry(round(KKm(:,1)*fs),1),'r*')
yline(ax2(2),Jm,'r','linewidth',1);
yline(ax2(2),-Jm,'r','linewidth',1);
yline(ax2(2),0,'k','linewidth',1);
set(ax2(1),'ylim',[0 max(SwimSp)]);
set(ax2(2),'ylim',[min(pry(Phase==-1,1)), max(pry(Phase==1,1))]);
legend('speed','body rotations')
linkaxes([sp1 sp2 ax2], 'x'); % links x axes

%HP acceleration - not enough memory for kk crossings, replaced with ylines
sp3 = subplot(6,1,6);
ax3= plotyy(k,SwimSp,k,Ahf(:,3));
set(ax3(2),'nextplot','add')
yline(ax3(2),Ja,'r','linewidth',1);
yline(ax3(2),-Ja,'r','linewidth',1);
yline(ax3(2),0,'k','linewidth',1);
set(ax3(1),'ylim',[0 max(SwimSp)]);
set(ax3(2),'ylim',[min(Ahf(Phase==-1,1)), max(Ahf(Phase==1,1))]);
legend('speed','Highpass Acceleration')
linkaxes([sp1 sp3 ax3], 'x'); % links x axes

% final diagnostic check, inspect the raw and smoothed pitch data to make
% sure low frequency/aplitude stroking hasn't been missed. 
figure(10); clf; 
sp1 = subplot(5,1,1:3);
plot(k,p,'k','Linewidth',3) ; hold on; set(gca,'ydir','rev','ylim',[min(p) max(p)]);
hold on
plot(t*fs,agree_yes,'-y','Linewidth',3);
plot(t*fs,onlym,'-m','Linewidth',3);
plot(t*fs,onlya,'-r','Linewidth',3);
legend('Track','Agree','Only Rot' ,'Only Acc');

sp2= subplot(5,1,4);
ax2= plotyy(k,SwimSp,k,pry(:,1));
set(ax2(2),'nextplot','add')
plot(ax2(2),KKm(:,1)*fs,pry(round(KKm(:,1)*fs),1),'r*')
yline(ax2(2),Jm,'r','linewidth',1);
yline(ax2(2),-Jm,'r','linewidth',1);
yline(ax2(2),0,'k','linewidth',1);
set(ax2(1),'ylim',[0 max(SwimSp)]);
set(ax2(2),'ylim',[min(pry(Phase==-1,1)), max(pry(Phase==1,1))]);
legend('speed','body rotations')
sp3 = subplot(5,1,5);
plot((1:length(p)),pitch*180/pi)
hold on
plot((1:length(p)),smoothpitch*180/pi,'r')
legend('pitch','smoothpitch')
linkaxes([sp1 sp2 sp3 ax2], 'x'); % links x axes


%% 8 MAKE 5SEC SUB-GLIDES
final_GDM = 'Mag'

if final_GDM == 'Acc'
    %SGtype indicates whether it is stroking (1) or gliding(0)
    glk1= ~glka; % AB - corrected
    SGtype=glk1;
    dur=5;
    SGL=splitGL_AB(dur,GLa);               % ES - changed from GL to GLa when using accel. method; AB - edited function: added floor() so remainders are ignored.
    SGLT=[SGL(:,1),SGL(:,2)-SGL(:,1)];
elseif final_GDM == 'Mag'
    %SGtype indicates whether it is stroking (1) or gliding(0)
    glk1= ~glkm; % AB - corrected
    SGtype=glk1;
    dur=5;
    SGL=splitGL_AB(dur,GLm);               % ES - changed from GL to GLa when using accel. method;  AB - edited function: added floor() so remainders are ignored.
    SGLT=[SGL(:,1),SGL(:,2)-SGL(:,1)];
end 

%check that all subglides have a duration of 5 seconds
% rangeSGLT=[min(SGLT(:,2)),max(SGLT(:,2))];
% delete(h3);
% figure(7);
% glk=eventon(SGLT,t);
% pgl=p;
% pgl(glk==0)=NaN;
% %h3=plot(SGL,t*fs,pgl,'m','Linewidth',3);


%% 9 Create summary table required for body density  

k=1:length(p);


%[MagAcc,pry,Sa,GL,KK] = magnet_rot_sa(Aw,Mw,fs,fluke_rate,f,alpha,1,k,Jm,tmax);
[smoothhead] = m2h(MagAcc.Mnlf(k,:),smoothpitch,smoothroll);
smoothpitchdeg=pitch*180/pi; %AB - use raw pitch
smoothheaddeg=smoothhead*180/pi;  % ES - added this line to convert from radians to degrees
smoothrolldeg=smoothroll*180/pi;  % ES - added this line to convert from radians to degrees

% ES - new Glide creation code, including:
% 1. Fix for acceleration values (i.e. in m/s2 not m/s/#samples/s) 
% 2. Fix to ensure acceleration calculated only if speed values available throughout glide (as in Miller 2016 Igor Pro procedures)
% 3. Fix for correct units of pitch, roll and heading
Glide = zeros(length(SGL),24) ;
for i=1:length(SGL)
    cue1=SGL(i,1)*fs;
    cue2=SGL(i,2)*fs;
    Glide(i,1)=SGL(i,1);%sub-Glide start point in seconds
    Glide(i,2)=SGL(i,2);%sub-Glide end point in seconds
    Glide(i,3)=SGL(i,2)-SGL(i,1);%sub-Glide duration
    Glide(i,4)=mean(p(round(cue1):round(cue2)));%mean depth(m)during sub-Glide
    Glide(i,5)=abs(p(round(cue1))-(p (round(cue2))));%total depth(m)change during sub-Glide
    Glide(i,6)=mean(SwimSp(round(cue1):round(cue2)));%mean swim speed during the sub-Glide, only given if pitch>30 degrees
    Glide(i,7)=mean((smoothpitchdeg(round(cue1):round(cue2))));%mean pitch during the sub-Glide
    Glide(i,8)=sin(mean(smoothpitch(round(cue1):round(cue2))));     % ES - changed this from "smoothpitchdeg" to "smoothpitch",i.e. changed from degrees to radians %mean sin pitch during the sub-Glide
    Glide(i,9)=std(smoothpitch(round(cue1):round(cue2)));           % ES - changed this from "smoothpitchdeg" to "smoothpitch",i.e. changed from degrees to radians  %SD of pitch during the sub-Glide
    %Glide(i,10)=mean(tempw(round(cue1):round(cue2)));               % ES - changed this from temp to tempw so refers to temperature around the whale not CTD temp %mean temperature during the sub-Glide
    %Glide(i,11)=mean(Dsw(round(cue1):round(cue2)));%mean seawater density (kg/m^3) during the sub-Glide
    
    if all(isnan(SwimSp(round(cue1):round(cue2))))==0;              % ES - added this line, if all speed values within the glide are not NaN
        if isnan(Glide(i,6)) ==0                                        % ES - added this line, so that only calcuate accelertion if there is a mean swim speed (otherwise accel given as 0 and following errors occur: "> In regress (line 84) Warning: X is rank deficient to within machine precision."
            try
                xpoly=(round(cue1):round(cue2))';
                % xpoly=(SGL(i,1):0.2:SGL(i,2))';
                ypoly=SwimSp(round(cue1):round(cue2));
                [B,BINT,R,RINT,STATS] = regress(ypoly,[xpoly ones(length(ypoly),1)]);
                
                Glide(i,12)=B(1)*fs;    % ES - changed this to ensure accel in m/s2, not m/fs %mean acceleration during the sub-Glide
                Glide(i,13)=STATS(1);   % ES - general note: this does not change with different gradients and SE %R2-value for the regression swim speed vs. time during the sub-Glide
                Glide(i,14)=STATS(4)*fs;% ES - changed this to correct for slope in s, in line with new gradient in ms2 %SE of the gradient for the regression swim speed vs. time during the sub-Glide;
            catch
                Glide(i,12)=NaN;%mean acceleration during the sub-Glide
                Glide(i,13)=NaN;%R2-value for the regression swim speed vs. time during the sub-Glide
                Glide(i,14)=NaN;%SE of the gradient for the regression swim speed vs. time during the sub-Glide
            end
        else    % ES - added these lines (to next end) to complete the if statement re. mean speeds not NaN
            Glide(i,12)=NaN;%mean acceleration during the sub-Glide
            Glide(i,13)=NaN;%R2-value for the regression swim speed vs. time during the sub-Glide
            Glide(i,14)=NaN;%SE of the gradient for the regression swim speed vs. time during the sub-Glide
        end
    else    % ES - added these lines (to next end) to complete the if statement re. SwimSp nan
        Glide(i,12)=NaN;%mean acceleration during the sub-Glide
        Glide(i,13)=NaN;%R2-value for the regression swim speed vs. time during the sub-Glide
        Glide(i,14)=NaN;%SE of the gradient for the regression swim speed vs. time during the sub-Glide
    end     % ES - added this to complete if statement
    sumphase=sum(Phase(round(cue1):round(cue2)));
    sp=NaN;
    sp(sumphase<0)=-1;
    sp(sumphase==0)= 0 ;
    sp(sumphase>0)=1 ;
    Glide(i,15)=sp;%Dive phase:0 bottom, -1 descent, 1 ascent, NaN not dive phase
    Dinf=D(find((D(:,1)*fs)<cue1 &(D(:,2)*fs)>cue2),:);
    if isempty(Dinf),
        Dinf=NaN(size(D,1),size(D,2)) ;
    end
    Glide(i,16)= Dinf(7);%Dive number in which the sub-Glide recorded
    Glide(i,17)=Dinf(6);%Maximum dive depth (m) of the dive
    Glide(i,18)=Dinf(3);%Dive duration (s) of the dive
    Glide(i,19)=circ_mean(smoothpitchdeg(round(cue1):round(cue2)));     %  ES - changed this from "smoothpitch" to "smoothpitchdeg" %Mean pitch(deg) calculated using circular statistics
    Glide(i,20)=1-(circ_var(smoothpitch(round(cue1):round(cue2))));         %Measure of concentration (r) of pitch during the sub-Glide (i.e. 0 for random direction, 1 for unidirectional)
    Glide(i,21)=circ_mean(smoothrolldeg(round(cue1):round(cue2)));      % ES - changed this from "smoothroll" to "smoothrolldeg" %Mean roll (deg) calculated using circular statistics
    Glide(i,22)=1-(circ_var(smoothroll(round(cue1):round(cue2))));          %Measure of concentration (r) of roll during the sub-Glide
    Glide(i,23)=circ_mean(smoothheaddeg(round(cue1):round(cue2)));      % ES - changed this from "smoothhead" to "smoothheaddeg" %Mean heading (deg) calculated using circular statistics
    Glide(i,24)=1-(circ_var(smoothhead(round(cue1):round(cue2))));          %Measure of concentration (r) of heading during the sub-Glide
    Glide(i,25) = SGL(i,3); % AB - subglide index

end
%% 10. Calculate glide ratio
% AB - Changed so it doesnt count glide or stroking time from excluded
% periods, addded max_depth
G_ratio = zeros(nn,14);            % ES - changed this from "G_ratio = zeros(size(D,1),10) ;" so as not to include last dive
for dive=  1:nnn                      % ES - changed this from "for dive=1:size(D,1)" so as not to include last dive
    kkdes=round(fs*D(dive,1)):round(fs*Bottom(dive,1)); %
    kkas=round(fs*Bottom(dive,3)):round(fs*D(dive,2)); % it is selecting the whole dive
    G_ratio(dive,1)=sum(Phase(kkdes)==-1)/fs;%total duration of the descet phase (s)
    G_ratio(dive,2)=length(find(SGtype(kkdes)==0))/fs;%total glide duration during the descent phase (s)
    G_ratio(dive,3)=G_ratio(dive,2)/G_ratio(dive,1);%glide ratio during the descent phase
    G_ratio(dive,6)=sum(Phase(kkas)==1)/fs;%total duration of the ascet phase (s)
    G_ratio(dive,7)=length(find(SGtype(kkas)==0))/fs;%total glide duration during the ascet phase (s)
    G_ratio(dive,8)=G_ratio(dive,7)/G_ratio(dive,6);%glide ratio during the ascent phase
    G_ratio(dive,11)=sum(Phase(kkdes)==0)/fs; %time excluded des
    G_ratio(dive,12)=sum(Phase(kkas)==0)/fs; %time excluded asc
    G_ratio(dive,13) = dive
    G_ratio(dive,14) = D(dive,6) %dive max depth
    if ~ ismember(dive,dives_w_exclusions) %AB - exclude from interrupted dive
        G_ratio(dive,4)=mean(smoothpitchdeg(kkdes));%mean pitch during the descet phase(degrees)     % ES - changed this from "smoothpitch" to "smoothpitchdeg"
        G_ratio(dive,5)=Bottom(dive,2)/G_ratio(dive,1);%descent rate (m/s)
        G_ratio(dive,9)=mean(smoothpitchdeg(kkas));%mean pitch during the ascent phase(degrees)       % ES - changed this from "smoothpitch" to "smoothpitchdeg"
        G_ratio(dive,10)=Bottom(dive,4)/G_ratio(dive,6);%ascent rate (m/s)
    end
end
%% Diagnostic Checks

% AB - Create the subset of data used in the MCMC analysis.
% Check expected relationship with depth. Look for noise and outliers. 
subs = Glide(:,[1,4,6,7,11,12,13,15,17,20,22]);
tab = array2table(subs, 'VariableNames',{'Start','Depth','Speed','Pitch','SwDens','Accel','R2','Phase','Max_Depth','pitch_conc','Roll_conc'});
log = abs(tab.Pitch) > 30 & tab.Max_Depth > 100 & abs(tab.Phase)>0 & ~isnan(tab.Accel) & ~isnan(tab.Speed) & ~isnan(tab.SwDens) & tab.Roll_conc > .9;
tab = tab(log,:);
c = tab.Pitch;
s = tab.Speed.^2*20; % marker size is proportional to the square of the speed and therefore the relative drag term.
figure
scatter(tab.Depth,tab.Accel,s,c,'filled')
title(join("n = "+height(tab)+final_GDM))
xlabel('Depth (m)') 
ylabel('Acceleration')
%ylim([-0.08,0.08])
legend('Speed^2')
c = colorbar;
c.Label.String = 'Pitch (degrees)';

%inspect individual glides in context
glides = 8%:length(tab.Accel)
figure
hold on
for n = glides
    hold on
    start = round(fs*tab.Start(n))-5*fs;
    glide = start:round(start+15*fs);
    t=(0:length(p(glide))-1)/fs;
    plot(t,(pry(glide,3)));
end
%ylim([-0.05,0.05])
yline(Ja,'r','linewidth',3)
xline(5)
xline(10)
yline(-Ja,'r','linewidth',3)


% 11. Export csvs
% ES - Changed from csvwrtie to dlmwrite. Do not use "csvwrite" as this reduces precision of results
filepath = fullfile(Glide_Sum_Dir,strcat(tag,'.csv'))
dlmwrite(filepath,Glide,'precision','%.7f')
filepath = fullfile(Glide_Ratio_Dir,strcat(tag,'.csv'))
dlmwrite(filepath,G_ratio,'precision','%.7f')
