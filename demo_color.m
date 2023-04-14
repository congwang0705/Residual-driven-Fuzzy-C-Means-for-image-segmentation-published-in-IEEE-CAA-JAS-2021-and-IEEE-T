%%
clc; clear all; close all;
%%
%Parameters
C = 2;  %cluster numbers
error = 1e-6;
m=2;
se=5;
%%
I = imread('86016.jpg');
figure(1);
imshow(I)

alpha=3;
If=w_ColorRecons_CO(double(I),se);
Isum=(double(I)+alpha*double(If))./(1+alpha);

II=colorspace('Lab<-RGB',uint8(Isum));
II = double(II);

[height,width,~] = size(I);
II1 = II(:,:,1);
II2 = II(:,:,2);
II3 = II(:,:,3);
lambda = zeros(3,1);
lambda(1) = 0*std(II1(:));%this parameter 1/4 could be tuned
lambda(2) = 0*std(II2(:));%this parameter 1/4 could be tuned
lambda(3) = 0*std(II3(:));%this parameter 1/4 could be tuned

para=4.6e1;
lambda1 = zeros(3,1);
lambda1(1)=para*std(II1(:));
lambda1(2)=para*std(II2(:));
lambda1(3)=para*std(II3(:));
lambda2 = zeros(3,1);
lambda2(1) = 1;%this parameter 1/4 could be tuned
lambda2(2) = 1;%this parameter 1/4 could be tuned
lambda2(3) = 1;%this parameter 1/4 could be tuned

tic;
[E,center,U] = RFCM_color(II,C,error,m,lambda,lambda1,lambda2);
toc;

[~, label] = max(U', [], 2);
label = reshape(label,height,width);
label=w_recons_CO(double(label),strel('square',se));
center_l=center(:,1);center_a=center(:,2);center_b=center(:,3);center_lab=cat(3,center_l,center_a,center_b);
centernew=255*colorspace('RGB<-Lab',center_lab);
%centernew=center_lab;
fs=reshape(centernew(label,:), height,width, 3);
figure(2);
imshow(uint8(fs),'border','tight')