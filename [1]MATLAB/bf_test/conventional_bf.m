%% Beamforming
clc;close all

%% Initialization
fc = 1e9;
lambda = physconst('LightSpeed')/fc;
array = phased.ULA('NumElements',10,'ElementSpacing',lambda/2);
array.Element.FrequencyRange = [8e8 1.2e9];

%% Rectangular pulse
close all

t = linspace(0,0.3,300)';
testing = zeros(size(t));
testing(201:205) = 1; %  A simple rectangular pulse.
% Pulse incidnet on the ULA from the following azumuth angle
ang_arr = [30; 0];
% The same signal impinges on 10 different antenna elements.
% With delay between adjacent antenna element, the signal is different at
% each antenna element.
x = collectPlaneWave(array,testing,ang_arr,fc);
% Add complex-valued Gaussian noise to the singal x
rng default
% npower is the variance. I have no clue why power is set to 0.5. 
npower = 0.5; 
x = x + sqrt(npower/2)*(randn(size(x)) + 1i*randn(size(x)));
% xind = x>0.8;
subplot(221)
plot(t,abs(x(:,1)))
title('Element 1 (magnitude)')
axis tight
ylabel('magnitude')
subplot(222)
plot(t,abs(x(:,2)))
title('Element 2 (magnitude)')
axis tight
ylabel('magnitude')
subplot(223)
plot(t,abs(x(:,3)))
title('Element 3 (magnitude)')
axis tight
ylabel('magnitude')
subplot(224)
plot(t,abs(x(:,4)))
title('Element 4 (magnitude)')
axis tight
ylabel('magnitude')

%% Construct a phase-shift beamformer
beamformer = phased.PhaseShiftBeamformer('SensorArray',array,...
    'OperatingFrequency',fc,'Direction',ang_arr,...
    'WeightsOutputPort',true);
[y,w] = beamformer(x); %% ???????????????
subplot(211)
plot(t,abs(testing))
axis tight
title('Original Signal')
ylabel('Magnitude')
subplot(212)
plot(t,abs(y))
axis tight
title('Received Signal with Beamforming')
ylabel('Magnitude')
xlabel('Seconds')

%% Plot the array normalized power response with and without beamforming weights.
azang = -180:30:180;
subplot(211)
% Pattern plots the array response pattern with the calculated weights. To
% me, it shows the power pattern for angle from [-180: 180] when the weight
% of the antenna element is fixed to the input weight.
pattern(array,fc,[-180:180],0,'CoordinateSystem','rectangular',...
    'Type','powerdb','PropagationSpeed',physconst('LightSpeed'))
set(gca,'xtick',azang);
title('Array Response without Beamforming Weights')
subplot(212)
pattern(array,fc,[-180:180],0,'CoordinateSystem','rectangular',...
    'Type','powerdb','PropagationSpeed',physconst('LightSpeed'),...
    'Weights',w)
set(gca,'xtick',azang);
title('Array Response with Beamforming Weights')