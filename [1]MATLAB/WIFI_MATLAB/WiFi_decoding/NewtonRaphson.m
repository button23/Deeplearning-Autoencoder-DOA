clear
close all
%%
syms x
f = @(x) exp(-1/2*x^2)+0.1*exp(-1/2*(x-10)^2);
df=@(x) - (exp(-(x - 10)^2/2)*(x - 10))/10 - x*exp(-x^2/2);
ddf = @(x) (exp(-(x - 10)^2/2)*(x - 10)^2)/10 - exp(-(x - 10)^2/2)/10 - exp(-x^2/2) + x^2*exp(-x^2/2);
t = 11;
tol = 10^(-5);
maxn = 1000;
for i =1 : maxn
    if abs(ddf(t)) < tol
        fprintf('Try a better initial guess\n');
        return;
    end
    t1=t-df(t)/ddf(t);
        if abs(t1-t) < tol
            fprintf('The root has converged at x = %.10f\n', t1);
            return;
        end
    t = t1;
    if i == maxn
        fprintf('The root has not converged!\n');
    end
end
result = t;
%%
d = -3:0.01:13;
yy = zeros(length(d), 1);
for j = 1 : length(d)
    yy(j) =f(d(j));
end
figure(1)
plot(d, yy)
