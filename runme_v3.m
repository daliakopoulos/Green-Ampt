% Author: Ioannis Daliakopoulos 
% Contact: daliakopoulos@hydromech.gr
% Date: 2016/12/01
% Revision: 2.0
% Copyright: Ioannis Daliakopoulos 2016

% Based on the fortran code of J.E. Parsons and R. Mu?oz-Carpena 
% (provided in: http://abe.ufl.edu/carpena/software/wingampt.shtml). 
% Same as their progam, the code submitted here follows 
% the method of Chu (1978), Mein and Larson (1971, 1973) 
% and Skaggs and Kaheel (1982) to calculate infiltration 
% for unsteady rainfall using the Green-Ampt equation. 

clear all; close all; clc; 

deltim = .05;   % calculation dt [h]
timoff = 0;     % time offset [h]

%soil variables
vsatk = .044;   % Ks(cm/h)
sav = 22.4;     % Sav(cm)
wcsat = .499;   % Sat.WC,
wcini = 0.25;   % Init.WC (cm^3/cm^3)
stmax = 0.5;   % maximum surface storage s(cm)
%==============

% rain_vector: [start hours, end hours, rainfall intensity cm/h]
rain_vector = [0.0  1.0 1.5%  
               1.0  4.5 0.1
               4.5  6.5 0.5];

[ntimes, times, ttp, ttpp, rint, prec, bf, fpp, f, stor, ro] = ...
green_ampt_v2(deltim, vsatk, sav, wcsat, wcini, stmax, rain_vector);

% 
% -----------------------------------------------------------------------------
%    Time    tp    tpp     R       P       F       fp       f       S       RO
%     h       h     h     cm/h     cm      cm     cm/h    cm/h     cm       cm
%  -----------------------------------------------------------------------------
% The variables are 
%  
% Ks = vertical saturated conductivity (cm/h) 
% Sav = average suction across the wetting front (cm) 
% wcs = water content, theta at saturation (cm^3/cm^3) 
% wci = initial water content, theta (cm^3/cm^3) at the start 
% M = wcs - wci (cm^3/cm^3) 
% ln = natural logarithm function 
% Fp = cumulative infiltration at ponding (cm) 
% fp = infiltration rate at ponding (cm/h) 
% tp = time to ponding (h) 
% tp’ = time required to infiltrate Fp if the system had started in ponded 
% conditions (h) 
% F = cumulative infiltration during the event (cm) 
% f = infiltration rate (cm/h) 
% t = time (h) 
% R = rainfall rate (cm/h) 
% P = amount of rainfall (cm) 
% S = surface storage (cm) 
% Smax = maximum surface storage (cm) 
% RO = runoff (cm) 
% dP, dS, dF = change in P, S, F since the last time step 

output = nan(ntimes, 10);

for ii=1:ntimes
    if ttp(ii) < 9999 || ttpp(ii) < 9999
       output(ii,:) = [times(ii),ttp(ii),ttpp(ii),rint(ii),prec(ii),bf(ii),fpp(ii),f(ii),stor(ii),ro(ii)];
    else
       output(ii,:) = [times(ii), NaN, NaN,rint(ii), prec(ii),bf(ii),fpp(ii),f(ii),stor(ii),ro(ii)];
    end
end


%%
figure;
subplot(3,1,1:2);
[ax, ~, ~] = plotyy(times,rint,times,stor);
axes(ax(1)); ylabel('Rainfall intensity [cm/h]');
axes(ax(2)); ylabel('Soil moisture [cm]');
xlabel('Time [h]');

subplot(3,1,3);
[ax, ~, ~] = plotyy(times,f,times,ro);
axes(ax(1));ylabel('Infiltration rate [cm/h]');
axes(ax(2));ylabel('Runoff [cm/h]');
xlabel('Time [h]');
axis tight;