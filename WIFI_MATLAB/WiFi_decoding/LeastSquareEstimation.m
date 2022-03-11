%% Least Squares Estimate
rand('state',100); % initializing the random number generation
y = [5:3:50]; % observations, y_i
y = y + 5*rand(size(y)); % y_i with noise added
x = 1:length(y); % the x co-ordinates

% Formulating in matrix for solving for least squares estimate
Y = y.';
X = [x.' ones(1,length(x)).'];
alpha = inv(X'*X)*X'*Y; % solving for m and c

% constructing the straight line using the estimated slope and constant
yEst = alpha(1)*x + alpha(2);

close all
figure
plot(x,y,'r.')
hold on
plot(x,yEst,'b')
legend('observations', 'estimated straight line')
grid on
ylabel('observations')
xlabel('x axis')
title('least squares straight line fit')

%% QR factorization
A = [1,1,1,1;-1,4,4,-1;4,-2,2,0]';
b = [1,2,3,4]';
[Q,R] = qr(A,0);
z = Q.'*b;
x_obs = inv(R)*z;

x_ver = real(R \z);

%%