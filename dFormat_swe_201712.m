%% dFormat_WL.m: Blabla
% Description
%%
clear all; close all; clc

%% Load and format tidal river water level data
% <html><h3>Load PVCA02</h3></html>
%
% PVCA02 is an object with raw surface water elevation (SWE) data measured 
% at PVCA02 station (table _d0_ ), and in which the processed data (table _dSWE_) 
% are stored
%
% * d0.tNum: time as a Matlab serial time
% * d0.tStr: stime in string format
% * d0.swe: surface water elevation w/r datum
% * d0.offset: sensor level w/r datum
% * d0.temp: water temperature

load Data/PVCA02

% Create a new table _d_ to store d0 and further processed data
d = d0;

% Water level d0.swe is in m w/r to datum
% Transform data w/r to river bed by adding 1 meter(Observation-based approximation)
d.swe_dat = d.swe; % swe w/r datum
d.swe = d.swe_dat+1; % swe w/r river bed 

% Start and end time of periods at which values are missing (failing sensor)
tFS1 = '12/07/2014  8:30:00 AM'; % Start time period 1
tFE1 = '12/09/2014  3:15:00 PM'; % End time period 1
tFS2 = '10/29/2015  8:30:00 AM'; % Start time period 2
tFE2 = '12/16/2015  3:15:00 PM'; % End time period 2

%%
% <html><h3>Find raw indices when values are missing</h3></html>
%
i_MV1 = find(d.tNum>= datenum(tFS1) & d.tNum<= datenum(tFE1)); 
i_MV2 = find(d.tNum>= datenum(tFS2) & d.tNum<= datenum(tFE2));
%%
% <html><h3>Cublic spline interpolation</h3></html>
%
% Fill swe data where values are missing because water level dropped under 
% the sensor level
d.swe_s3 = spline(d.tNum, d.swe, d.tNum);

% Replace by NaN data missed by interpolation 
d.swe_s3([i_MV1; i_MV2]) = NaN;
d.swe_s3(end-15 : end) = NaN;

% Display data before and after interpolation
figure
subplot 211 % all the data
plot(d.tNum, d.swe_s3,'r.'); hold on
plot(d.tNum, d.swe,'b.')
datetick('x', 12)
subplot 212 % select a data window to have a closer look at the reconstruted data
plot(d.tNum(1:500), d.swe_s3(1:500),'r'); hold on
plot(d.tNum(1:500), d.swe(1:500),'b')
legend('raw','intrapolated')

%% Fill missing values using harmonic decomposition and reconstruction
% Filling missing values(_MV1_ and _MV2_) is a necessary step to avoid edge 
% effects during the complex demodulation. We do so by using harmonic decomposition 
% with the _UTide_ Matlab functions (Codiga, 2011). Theis package needs to be 
% saved in workgin directory, or pointed to in order for this script to run.
% The reconstructed values will be replaced by NaNs after complex demodulation 
% is executed.
% 
% Codiga, D.L., 2011. Unified Tidal Analysis and Prediction Using the UTide 
% Matlab Functions. Technical Report 2011-01. Graduate School of Oceanography, 
% University of Rhode Island, Narragansett, RI. 59pp. 
% ftp://www.po.gso.uri.edu/pub/downloads/codiga/pubs/2011Codiga-UTide-Report.pdf

% Transfrom time to GMT time as required by _UTide_ 
tGMT = d.tNum + 6*1/24; 
swe_M = d.swe_s3; % Measured SWE

% Select one year of data as training data
tS4t = datenum('5/15/2013'); 
tE4t = datenum('5/24/2014');
i_d4t = find(tGMT >=tS4t & tGMT<=tE4t); 
t4t = tGMT(i_d4t); % time for training data
swe4t = swe_M(i_d4t); % swe data for training

% Harmonic decomposition on the training data with _ut_solve_ function
coef = ut_solv(t4t,swe4t,[],10, 'auto','White','OLS','LinCI','DiagnPlots',...
    'OrderCnstit' , 'frq','NoTrend');  

% Reconstruction over the period corresponding to MVs with _ut_reconstr_ function
[swe_R1, ~ ] = ut_reconstr(tGMT(i_MV1), coef);
[swe_R2, ~ ] = ut_reconstr(tGMT(i_MV2), coef);

d.sweF = d.swe; % SWE Filled
d.swe_s3F = d.swe_s3; % SWE intrapolated and filled
d.sweF(i_MV1) = swe_R1; 
d.sweF(i_MV2) = swe_R2;
d.swe_s3F(i_MV1) = swe_R1;
d.swe_s3F(i_MV2) = swe_R2;

% Display data before and after filling
figure; hold on
plot(tGMT, d.swe_s3F,'r')
plot(tGMT, d.swe_s3,'b')
legend('Initial data','Reconsructed data')

%% Complex demodulation and low-pass filter
% Extract temporal variation of amplitude of major harmonic constituents
t = d.tNum;
swe = d.swe_s3F;
si = 15; % data sampling interval in minutes

%%
% <html><h3>4-hour demodulation</h3></html>
%
demodPer = 4; % demodulation period in hours
plpf = 8; % period of low pass filter in hours
swe4 = compDem(demodPer, t, swe, plpf, si); % amplitude of 4h-demodulated signal 
[swe4lp,coef] = lanczos(swe4 , 24*30, si); % amplitude of 4h-demodulated signal after low-pass filter

% Replace filtered data corresponding to the MVs by NaNs
swe4([i_MV1' i_MV2']) = NaN; 
swe4lp([i_MV1' i_MV2']) = NaN;

% Results Visualization
figCmpMod(t, swe, swe4, swe4lp, 4)
%%
% <html><h3>6-hour demodulation (Quarter-diurnal)</h3></html>
%
demodPer = 6;
plpf = 13;
swe6 = compDem(demodPer, t, swe, plpf, si);
[swe6lp,coef] =  lanczos(swe6 , 24*30, si);

% Replace filtered data corresponding to the MVs 
swe6([i_MV1' i_MV2']) = NaN; 
swe6lp([i_MV1' i_MV2']) = NaN;

% Results Visualization
figCmpMod(t, swe, swe6, swe6lp, 6)

%%
% <html><h3>12-hour demodulation (Semi-diurnal)</h3></html>
%
demodPer = 12;
plpf = demodPer *2;
swe12 = compDem(demodPer, t, swe, plpf, si);
[swe12lp,~] =  lanczos(swe12, 24*30, si);

% Replace filtered data corresponding to the MVs 
swe12([i_MV1' i_MV2']) = NaN; 
swe12lp([i_MV1' i_MV2']) = NaN;

% Results Visualization
figCmpMod(t, swe, swe12, swe12lp, 12)

%%
% <html><h3>24-hour demodulation (Diurnal)</h3></html>
%
demodPer = 24;
plpf = 50;
swe24 = compDem(demodPer, t, swe , plpf, si);
[swe24lp,~] =  lanczos(swe24 , 24*30, si);

% Replace filtered data corresponding to the MVs 
swe24([i_MV1' i_MV2(1)-110:i_MV2(end)+110]) = NaN; 
swe24lp([i_MV1' i_MV2']) = NaN;

% Results Visualization
figCmpMod(t, swe, swe24, swe24lp, 24)

%%
% <html><h3>30-hour demodulation </h3></html>
%
demodPer = 30;
plpf =60;
swe30 = compDem(demodPer, t, swe , plpf, si);
[swe30lp,~] =  lanczos(swe30 , 24*15, si);

% Replace filtered data corresponding to the MVs 
swe30([i_MV1' i_MV2(1)-110:i_MV2(end)+110]) = NaN; 
swe30lp([i_MV1' i_MV2']) = NaN;

% Results Visualization
figCmpMod(t, swe, swe30, swe30lp, 30)

%%
% <html><h3>13.6625-day demodulation (Fortnightly)</h3></html>
%
demodPer = 24*13.66;
plpf = 24*30;
swe15j = compDem(demodPer, t, swe, plpf, si);
% Replace filtered data corresponding to the MVs 
swe15j([i_MV1; i_MV2]) = NaN;

figure; hold on
plot(t, swe, 'k')
plot(t, swe15j, 'k')

%%
% <html><h3>27.5 days demodulation (MM harmonic constituent)</h3></html>
%
demodPer = 661.3;% 27.5j
plpf = 24*60;
sweMM = compDem(demodPer,t,swe , plpf, si);

% Replace filtered data corresponding to the MVs 
sweMM([i_MV1; i_MV2]) = NaN;

figure; hold on
plot(t, sweMM, 'k')
close

% Subtidal as defined in Henrie et al., 2014
wl30hLlp = lanczos(swe-nanmean(swe), 30,si);

figure; hold on
plot(t, d.wls3)
plot(t, wl30hLlp, 'k')

close

% Subtidal: 
% 1. Filter out frequencies lower than 1/30day. NOT complex
% demodulation
swe15Llp = lanczos(swe-nanmean(swe), 24*15,si);
swe30Llp = lanczos(swe-nanmean(swe), 24*30,si);
sweSubtid15Llp([i_MV1; i_MV2]) = NaN;
sweSubtid30Llp([i_MV1; i_MV2]) = NaN;

% 2. Reconstruct signal with SSA and SA constituents
tGMT = t + 6*1/24;  % GMT time required by U_solve (GMT time is CR time +6 hours )
swe_M = d.wls3-nanmean(d.wls3);
lat = 10;
cstt = {'SSA', 'SA'}; 
coef = ut_solv(tGMT,swe_M,[],lat, cstt,'White','OLS','LinCI','DiagnPlots',...
    'OrderCnstit' , 'frq','NoTrend');  
sweSASSA = ut_reconstr(tGMT, coef);
 
% 3. Substract SASSA from 30Llp
sweSubtid15 = swe15Llp - sweSASSA;
sweSubtid15lp = lanczos(sweSubtid15, 24*30, 15);
sweSubtid15lp([i_MV1; i_MV2]) = NaN;

sweSubtid30 = swe30Llp - sweSASSA;
sweSubtid30([i_MV1; i_MV2]) = NaN;

figure; hold on
plot(t, swe15Llp, 'k')
plot(t, swe30Llp, 'k')
plot(t, sweSASSA, 'b')
plot(t, sweSubtid15, 'r')
plot(t, sweSubtid30, 'r')

figure; hold on
plot(swe_M)
plot(sweSASSA)
plot([1 111438], [nanmean(swe_M), nanmean(swe_M)])
% close

% Save all constituents as additional column in the table
d.swe4 = swe4; 
d.swe4lp = swe4lp;
d.swe6 = swe6; 
d.swe6lp = swe6lp;
d.swe12 = swe12; 
d.swe12lp = swe12lp;
d.swe24 = swe24; 
d.swe24lp = swe24lp;
d.swe15j = swe15j;
d.sweMM = sweMM;

d.swe15Llp = swe15Llp;
d.swe30Llp = swe30Llp;
d.sweSASSA = sweSASSA;
d.sweSubtid15 = sweSubtid15;
d.sweSubtid15lp = sweSubtid15lp;

% Sum of small frequencies harmonics amplitudes
d.sumSfq = d.swe4lp + d.swe6lp + d.swe12lp + d.swe24lp;

swe = d;

save ~/Dropbox/WORK/TBProject/Data/PVNP/frMat/PVCA02_Corr d -append
save ~/Dropbox/WORK/TBProject/analysis/TidalInf/Final/Data/dwl swe i_MV1 i_MV2 -append
%% --- DATA STANDARDIZATION: Unity-based normalization
swe = d;
sweN = table; 
sweN.tNum = swe.tNum; sweN.tStr = swe.tStr;

varN = swe.Properties.VariableNames;
for i = 11:18
    var2N = eval(['wl.',varN{i}]); 
    vN =  (var2N - min(var2N))/(max(var2N) - min(var2N));
    eval(sprintf('%s = vN', ['wlN.', varN{i}]));
end

save Data/dwl sweN -append

