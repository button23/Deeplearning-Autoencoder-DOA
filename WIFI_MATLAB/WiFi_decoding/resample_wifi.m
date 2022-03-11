clc;
close all;

% ff = input('Enter Frequency')
% ffs = input('Enter Sampling Frequency')

ff = 100;
ffs = 200;

% Specify Time Range of Signal
t = 0:0.1:100
ts = 0:10:100;

x = sin(2*3.14*ff*t);
subplot(2,1,1)
plot(t,x)

% Write sampling formula
y = sin(2*3.14*ff*ts/ffs);
subplot(2,1,2);
stem(ts,y)