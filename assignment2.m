%% 1
clear all;
close all;

f = imread("barbara_512.png");
[m, n] = size(f);
a = zeros(m, n);
b = zeros(m, n);
c = zeros(m, n);
d = zeros(m, n);
e = zeros(m, n);

for i=2:m-1
    for j=2:n-1
        a(i, j) = f(i+1, j) - f(i, j);
        b(i, j) = f(i, j) - f(i-1, j);
        c(i, j) = (f(i+1, j)- f(i-1, j))/2;
        d(i, j) = f(i+1, j)- (2*f(i, j)) + f(i-1, j);
        e(i, j) = f(i, j+1)- (2*f(i, j)) + f(i, j-1);
    end
end

figure; imshow(f, []); title("original");
figure; imshow(a, []); title("a");
figure; imshow(b, []); title("b");
figure; imshow(c, []); title("c");
figure; imshow(d, []); title("d");
figure; imshow(e, []); title("e");

%% 2
clear all;
close all;

f = double(imread("barbara_512.png"));
[m, n] = size(f);
u = 0;
sigma = 1;
a = zeros(m, n);
b = zeros(m, n);
c = zeros(m, n);

L1=zeros(3,3);
L1(1, 1)=0; L1(1, 2)=1;  L1(1,3)=0;
L1(2, 1)=1; L1(2, 2)=-4;  L1(2, 3)=1; 
L1(3, 1)=0; L1(3, 2)=1;  L1(3, 3)=0;

L2=zeros(3,3);
L2(1, 1)=1; L2(1, 2)=0;  L2(1,3)=1;
L2(2, 1)=1; L2(2, 2)=-4;  L2(2, 3)=0; 
L2(3, 1)=1; L2(3, 2)=0;  L2(3, 3)=1;

Linf=zeros(3,3);
Linf(1, 1)=1; Linf(1, 2)=1;  Linf(1,3)=1;
Linf(2, 1)=1; Linf(2, 2)=-8;  Linf(2, 3)=1; 
Linf(3, 1)=1; Linf(3, 2)=1;  Linf(3, 3)=1;

for i=2:m-1
    for j=2:n-1
        a(i, j)=(L2(1, 1)*f(i-1, j-1)+L2(1, 2)*f(i-1, j)+L2(1,3)*f(i-1, j+1) ...
                +L2(2, 1)*f(i, j-1)+L2(2, 2)*f(i, j)+L2(2, 3)*f(i, j+1)...
                +L2(3, 1)*f(i+1, j-1)+L2(3, 2)*f(i+1, j)+L2(3, 3)*f(i+1, j+1));
        b(i, j)=(1/2)*(L1(1, 1)*f(i-1, j-1)+L1(1, 2)*f(i-1, j)+L1(1,3)*f(i-1, j+1) ...
                +L1(2, 1)*f(i, j-1)+L1(2, 2)*f(i, j)+L1(2, 3)*f(i, j+1)...
                +L1(3, 1)*f(i+1, j-1)+L1(3, 2)*f(i+1, j)+L1(3, 3)*f(i+1, j+1));
        c(i, j)=(Linf(1, 1)*f(i-1, j-1)+Linf(1, 2)*f(i-1, j)+Linf(1,3)*f(i-1, j+1) ...
                +Linf(2, 1)*f(i, j-1)+Linf(2, 2)*f(i, j)+Linf(2, 3)*f(i, j+1)...
                +Linf(3, 1)*f(i+1, j-1)+Linf(3, 2)*f(i+1, j)+Linf(3, 3)*f(i+1, j+1));
    end
end
figure; imshow(f, []); title("original");
figure; imshow(a, []); title("l 2");
figure; imshow(b, []); title("l 1");
figure; imshow(c, []); title("l infinity");


%% 3
clear all;
close all;

f = imread("noisy_pollens.png");
[m, n] = size(f);
a = zeros(m, n);

figure;
plot = histogram(f);
figure;
imshow(f, []);
%salt and pepper noise
figure; imshow(medfilt2(f), [0, 255]); title("denoised");
f = medfilt2(f);
%contrast scaling
g2=0.3;
u2 = uint8( 255 * (double(f)/255).^(1/g2) );

figure; imshow(u2, [0, 255]); title("with gamma correction");
