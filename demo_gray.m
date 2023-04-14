%%
clc; clear all; close all;
Io = imread('sy2.bmp');
rng('default');
[m, n, p] = size(Io);
if p~=1
    Io=rgb2gray(Io);
    p=1;
end
figure(1)
imshow(Io,'border','tight');
%%
I=imnoise(uint8(Io),'gaussian',0,15^2/255^2);
I=imnoise(uint8(I),'salt & pepper',0.15);
I=uint8(I);
figure(2)
imshow(I,'border','tight');
%%
[counts,x]=imhist(I);
figure(3)
imhist(I);
%%
se=3;
If=w_recons_CO(double(I),strel('square',se));
If=uint8(If);
figure(4)
imshow(If,'border','tight');
%%
alpha=0;
Isum=(double(I)+alpha*double(If))./(1+alpha);
Isum=uint8(Isum);
figure(5)
imshow(Isum,'border','tight');
%%
X = reshape(double(Isum), m*n, p);
%%
%Parameters
C = 4;error = 1e-6;b=2;
lambda = 0*std(X);
lambda1 = 7.45e-2*std(X); 
lambda2 = 0.0008; 

tic;
Noi=double(I-Io);Noi=Noi(:);
[E,W,center,U,iter,obj,J,diff] = RFCM(double(Isum),Noi,C,error,b,lambda,lambda1,lambda2);
toc;
[~,label]=max(U',[],2);
label=reshape(label, m, n, p);
label=w_recons_CO(double(label),strel('square',2));
Is=reshape(center(label(:), :), m, n, p);
Is=w_recons_CO(double(Is),strel('square',se));
figure(6)
imshow(uint8(Is),'border','tight')