%% Reinitialize
%clear
%clc
%close all

%% Import
load('minDt.mat');
load('EminDt.mat');
load('time.mat');
load('distance.mat');

%% Plot
% MinDt Plot
figure(8)
subplot(2,1,1)
surf(time,distance,minDt');
axis([0 max(time) 0 max(distance)]);
set(gca, 'FontSize',12);
xlabel('time(s)', 'FontSize',12);
ylabel('distance along major axis({\mu}m)','FontSize',12);
title('Density of minDt on Membrane (Ellipsoid, Original)','FontSize',16);
c=colorbar;
ylabel(c,'N_{minDt}(t)');

% EMinDt Plot
subplot(2,1,2)
surf(time,distance,EminDt');
axis([0 max(time) 0 max(distance)]);
set(gca, 'FontSize',12);
xlabel('time(s)', 'FontSize',12);
ylabel('distance along major axis({\mu}m)','FontSize',12);
title('Density of EminDt on Membrane (Ellipsoid, Original)','FontSize',16);
c=colorbar;
ylabel(c,'N_{EminDt}(t)');