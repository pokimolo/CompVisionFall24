clc; close all; clear all;

input = double(imread("/MATLAB Drive/Computer Vision/woman_blonde.png"));
SUP=max(max(input));
inputdivmax=input/SUP;
f = imnoise(inputdivmax,"gaussian",0, 0.002);
figure; imshow(f); title("noisy image f");

%Perona Malik with auto stop:
u = f; 
k = 0.02;
dt = 0.25;
[m, n] = size(f);
stop = estimate_var(f);

while true
    u = expand(expand(u)); 
    u_new = zeros(size(u));
    s = (Dx_0_circ(u)).^2 + (Dy_0_circ(u)).^2;
    alpha = 1 ./ (1 + (s / k^2));  
    
    u_new = u + dt * (Dx_minus_circ(alpha .* Dx_plus_circ(u)) + Dy_minus_circ(alpha .* Dy_plus_circ(u)));  
    u = shrink(shrink(u_new));  
    % (e) comparing up to 4 decimal places
    % Display the mean of the image u
    disp(['Mean of u: ', num2str(mean(u(:)))]);  % Mean of u

    % Display the mean of the image f
    disp(['Mean of f: ', num2str(mean(f(:)))]);  % Mean of f

    temp = sum((u - f).^2, 'all') / (m * n - 1); 
    
    if temp > stop
        break;
    end
end

figure;
imshow(u, []);
title('Perona Malik Image');


function U=expand(u, n)
if nargin==1
    n=1;
end
for i=1:n
    [M, N]=size(u);
    U=u;
    V=[U(1, :); U ; U(M, :)];
    U=[V(:, 1), V , V(:, N)];
    u=U;
end
end

function u=shrink(U, n)
if nargin==1
    n=1;
end
[M, N]=size(U);
u=U(1+n:M-n, 1+n:N-n);
end

function d=Dx_plus_circ(u)
u=expand(u);
d=circshift(u, -1, 1)-u;
d=shrink(d);
end
function d=Dx_minus_circ(u)
u=expand(u);
d= u-circshift(u, 1, 1);
d=shrink(d);
end

function d=Dy_minus_circ(u)
u=expand(u);
d=u-circshift(u, 1, 2);
d=shrink(d);
end

function d=Dx_0_circ(u)
u=expand(u);
d=(circshift(u, -1, 1)-circshift(u, 1, 1))/2;
d=shrink(d);
end

function d=Dy_plus_circ(u)
u=expand(u);
d=circshift(u, -1, 2)-u;
d=shrink(d);
end

function d=Dy_0_circ(u)
u=expand(u);
d=(circshift(u, -1, 2)-circshift(u, 1, 2))/2;
d=shrink(d);
end

% (f) estimate_var function
function var_est = estimate_var(f)
    w = [1, -2, 1;
         -2, 4, -2;
          1, -2, 1];
    [M, N] = size(f);
    u = f;
    for i = 2:M-1
        for j = 2:N-1
            u(i, j) = sum(sum(w .* f(i-1:i+1, j-1:j+1)));
        end
    end
    var_est = sum(sum(u.^2)) * (1 / (36 * (M - 2) * (N - 2)));
end
