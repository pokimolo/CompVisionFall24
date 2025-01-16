clc;
clear all;
close all;
f = double(imread("barbara_512.png")); %gets the pixles from the image
[M, N] = size(f);
figure; imshow(f,[]);
title("Origional Image")

%% 1(a). Scale the image so that the intensity lies in the range 0 to 1.

f_scaled = (f - min(f(:))) / (max(f(:)) - min(f(:)));

figure; imshow(f_scaled);
title("intensity lies in the range 0 to 1");

%% 1(b) Write functions expand and shrink that expands and shrinks the image f by one-pixel width.
function U=expand(u)
[M, N]=size(u);
U=u;
V=[U(1, :); U ; U(M, :)];
U=[V(:, 1), V , V(:, N)];
end  

function u=shrink(U)
[M, N]=size(U);
u=U(2:M-1, 2:N-1);
end

figure;imshow(expand(f_scaled));title("expanded");
figure;imshow(shrink(f_scaled));title("shrunk");

%% 1(d).
max_iterations=1;
dt=0.05;
u=f_scaled;

for n=1:max_iterations
    u=expand(u);
    u_x=u;
    for i=2:M-1
        for j=2:N-1
            u_x(i, j) = u(i, j) + dt * (u(i-1, j) + u(i+1, j) + ...
                u(i, j-1) + u(i, j+1));
            u_x(i, j) = u_x(i, j)/(1 + 4*dt);
        end
    end
    u=shrink(u_x);
    u_x=u;
end

figure();
imshow(u_x);
title("Q1 part D")
%% 1(f). & 1(g).
max_iterations=401;
dt=0.05;
u=f_scaled;

for n=1:max_iterations
    u=expand(u);
    u_x=u;
    for i=2:M-1
        for j=2:N-1
            u_x(i, j) = u(i, j) + dt * (u(i-1, j) + u(i+1, j) + ...
                u(i, j-1) + u(i, j+1));
            u_x(i, j) = u_x(i, j)/(1 + 4*dt);
        end
    end
    u=shrink(u_x);
    u_x=u;
    if(n==50)
        t=max(max(u-f)).^2;
        fprintf("Q1 part F at 50:\t"+t+"\n");
        t2=max(max(u-f));
        fprintf("Q1 part G at 50:\t"+t2+"\n");
    elseif(n==100)
        t=max(max(u-f)).^2;
        fprintf("Q1 part F at 100:\t"+t+"\n");
        t2=max(max(u-f));
        fprintf("Q1 part G at 100:\t"+t2+"\n");
    elseif(n==200)
        t=max(max(u-f)).^2;
        fprintf("Q1 part F at 200:\t"+t+"\n");
        t2=max(max(u-f));
        fprintf("Q1 part G at 200:\t"+t2+"\n");
    elseif(n==300)
        t=max(max(u-f)).^2;
        fprintf("Q1 part F at 300:\t"+t+"\n");
        t2=max(max(u-f));
        fprintf("Q1 part G at 300:\t"+t2+"\n");
    elseif(n==400)
        t=max(max(u-f)).^2;
        fprintf("Q1 part F at 400:\t"+t+"\n");
        t2=max(max(u-f));
        fprintf("Q1 part G at 400:\t"+t2+"\n");
    end
end