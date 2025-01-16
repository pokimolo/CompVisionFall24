clc;
clear all;
close all;

% 1
f = double(imread("barbara_512.png"));
figure; imshow(f,[]);
[M, N] = size(f);

% 2
figure;
plot(f(100,:));
title("100th row");
figure;
plot(f(:,100));
title("100th col");

%3 average value
sum= 0;
num_pixels = M*N;
for i=1:M
    for j=1:N
        sum = sum +f(i,j);
    end
end
average_value = sum/num_pixels;
% about 112

% 4 Apply mean filter
u = zeros(M,N);
num_pixels = M*N;
for i=2:M-1
    for j=2:N-1
        mean=(f(i-1, j-1)+f(i,j-1)+f(i+1,j-1)+ ...
            f(i-1,j)+f(i,j)+f(i+1,j)+ ...
            f(i-1, j+1)+f(i, j+1)+f(i+1, j+1))/9;
        u(i,j) = mean;
    end
end

figure;
imshow(f, []);
title('Original Image');

figure;
imshow(u, []);
title('Blurred Image');


%5 
figure;
v = f - u;
imshow(v, []);
title('V = F - U');
%6
w = [0 1 0;
    1 -4 1;
    0 1 0];


h=f;
for i=2:M-1
    for j=2:N-1
        h(i, j)=(w(1, 1)*f(i-1, j-1)+w(1, 2)*f(i-1, j)+w(1,3)*f(i-1, j+1) ...
                +w(2, 1)*f(i, j-1)+w(2, 2)*f(i, j)+w(2, 3)*f(i, j+1)...
                +w(3, 1)*f(i+1, j-1)+w(3, 2)*f(i+1, j)+w(3, 3)*f(i+1, j+1));
    end
end

diff=f-h;
figure; imshow(h,[]);
title("Convolution");
figure; imshow(diff,[]);
title("Difference");
%%

%7 Convolve
Px = [-1 -1 -1;
    0 0 0;
    1 1 1];
Py = [-1 0 1;
    -1 0 1;
    -1 0 1];

px=f;
py=f;
for i=2:M-1
    for j=2:N-1
        px(i, j)=(Px(1, 1)*f(i-1, j-1)+Px(1, 2)*f(i-1, j)+Px(1,3)*f(i-1, j+1) ...
                +Px(2, 1)*f(i, j-1)+Px(2, 2)*f(i, j)+Px(2, 3)*f(i, j+1)...
                +Px(3, 1)*f(i+1, j-1)+Px(3, 2)*f(i+1, j)+Px(3, 3)*f(i+1, j+1));
    end
end
figure; imshow(px,[]);
title("Convolution Px");

for i=2:M-1
    for j=2:N-1
        py(i, j)=(Py(1, 1)*f(i-1, j-1)+Py(1, 2)*f(i-1, j)+Py(1,3)*f(i-1, j+1) ...
                +Py(2, 1)*f(i, j-1)+Py(2, 2)*f(i, j)+Py(2, 3)*f(i, j+1)...
                +Py(3, 1)*f(i+1, j-1)+Py(3, 2)*f(i+1, j)+Py(3, 3)*f(i+1, j+1));
    end
end

figure; imshow(py,[]);
title("Convolution Py");

%% 

%8
A = [1 2 3;
    4 5 6];

[m, n]=size(A);
[U, S, V]=svd(A);

for i=1:m
    s(i)= S(i,i);
end
Vt=V.';
X=U*S*Vt;
% Checked that SVD worked

% S has a silent column for the third column, U only has 2 columns unlike S
% and V

% 9 Reduced SVD

[rU, rS, rV] = svd(A, 'econ');
rVt=rV.';
rX=rU*rS*rVt;
% Checked that rSVD worked
% rV is a different size from rU and rS because it's the only nonsquare matrix



% 10

[fm, fn]=size(f);
[fU, fS, fV]=svd(f);

for i=1:fm
    fs(i)= fS(i,i);
end
%%
R=200; % How many columns do you want to transfer?

newU=fU(1:fm, 1:R);
newS=double(zeros(R,R));
newV=fV(1:fn, 1:R);

for i=1:R
    newfS(i, i)=fs(i);
end

newf=newU*newfS*newV';
%%

figure;
imshow(newf,[]);
title("Compressed");

percentage_of_data_saved=(fm*fn-(fm*R+fn*R))*100/(fm*fn);

sprintf("Percentage of data saved is %3.2f", percentage_of_data_saved)
% 21.88% saved data