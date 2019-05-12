d  = sin(0.05*pi*[1:200] + 2*pi*rand);
figure;
plot(d);
g = randn(1, 200);
v1 = filter(1, [1 -0.8], g);
v2 = filter(1, [1 0.6], g);
x = d + v1;
figure;
plot(x);

p = 12;
Rv2 = covar(v2, p);
rxv2 = convm(x, p)' * convm(v2, p) / (length(x)-1);
W = rxv2(1, : ) / (Rv2);
v1hat = filter(W, 1, v2);
dhat = x - v1hat;
figure;
plot(dhat);

function X = convm(x, p)
N = length(x) + 2*p-2;
x = x(:);
xpad = [zeros(p-1, 1);x;zeros(p-1, 1)];
for i=1:p
    X(:, i) = xpad(p-i+1:N-i+1);
end
end

function R = covar(x, p)
X = x(:);
m = length(x);
X = x - ones(m, 1)*(sum(x)/m);
R = convm(x, p)' * convm(x, p)/(m-1);
end
