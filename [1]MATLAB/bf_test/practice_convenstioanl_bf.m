close all
clear,clc
c = physconst('LightSpeed');
fc = 2.45e9;
txarray = phased.ULA('NumElements',10,'ElementSpacing',c/fc/2);
txarray.Element.FrequencyRange = [2e9 3e9];

t = linspace(0,1,1000);
tx1 = sin(2*pi*5*t)';
tx2 = sin(2*pi*10*t)';
tx = [tx1 tx2];

angle = [20 40;0 0];

figure
subplot(211)
plot(t,tx1,'r')
title('Sinewave at frequency 5Hz')
ylabel('Amplitude')
subplot(212)
plot(t,tx2,'b')
title('Sinewave at frequency 10Hz')
ylabel('Amplitude')
xlabel('time (second)')

% plot(t,tx)
%%
y = collectPlaneWave(txarray,tx,angle,fc);
% figure
% subplot(221)
% plot(t,y(:,1));
% subplot(222)
% plot(t,y(:,2));
% subplot(223)
% plot(t,y(:,3));
% subplot(224)
% plot(t,y(:,4));


%%
close all
beamformer = phased.PhaseShiftBeamformer('SensorArray',txarray,...
    'Direction',angle,'OperatingFrequency',fc,...
    'WeightsOutputPort',true);
[yy,w] = beamformer(y);

subplot(311)
plot(t,tx1)
axis tight
title('Original Signal')
ylabel('Magnitude')
hold on 
plot(t,tx2)
legend('5Hz','10Hz')

subplot(312)
plot(t,abs(y(:,5)))
axis tight
title('Received Signal without Beamforming')
ylabel('Magnitude')
xlabel('Seconds')

subplot(313)
plot(t,abs(yy(:,1)))
axis tight
title('Received Signal with Beamforming')
ylabel('Magnitude')
xlabel('Seconds')
hold on
plot(t,abs(yy(:,2)))
legend('5Hz','10Hz')

%%
% subplot(211)
% pattern(txarray,fc,[-180:180],0,'CoordinateSystem','rectangular',...
%     'Type','powerdb','PropagationSpeed',physconst('LightSpeed'))
% set(gca,'xtick',-180:30:180)
% subplot(212)
% pattern(txarray,fc,[-180:180],0,'CoordinateSystem','rectangular',...
%     'Type','powerdb','PropagationSpeed',physconst('LightSpeed'),...
%     'Weights',w)
% set(gca,'xtick',-180:30:180)
% hold on
% xline(angle(1,:),'r','LineWidth',1)





