clc; 
clear all;
close all;
lambda = 1*10^-5;
img = im2double(imread("/Users/anthony/Pictures/images/image5.jpg"));
cluster = 2;
%figure()
% imshow(img);
img = img(:,:,1);
[h,w,c] = size(img);
pixel = reshape(img,[h*w,c]);
mean = mean(pixel);
std = std(pixel);
pixel = (pixel-mean)./std;
pixel = (pixel-mean);
[idx, class] = kmeans(pixel, cluster);
kmeansAns = reshape(idx,[h,w]);
%figure()
% imshow(kmeansAns/3)
%% initialize
int_cov = {};
initial_priors = zeros(cluster,1);
for i = 1:cluster
    data = [];
    for j=1:size(idx)
        if idx(j) == i
            data(end+1,:)=pixel(j,:);
        end
    end
    int_cov{end+1} = cov(data);%+lambda*diag(ones(1, c));
    initial_priors(i) = length(data)/length(idx);
end
save('data.mat')
load('data.mat')
%% EM iteration
priors = initial_priors;
u = class
Cov=int_cov
old_log_likelihood = -inf;

% priors = [0.3;0.3;0.4];
% u = zeros(3);
% cov = {};
% for i =1:cluster
%     cov{end+1} = cov(u);
% end

% cov = cov + lambda*diag(ones(1, length(cov)));
for j = 1:5000
    p = [];
    for i =1:cluster
        p(end+1,:) = priors(i)*mvnpdf(pixel(:,:),u(i,:),Cov{i});
    end
    p=p';
    p_x = sum(p,2);
    log_likelihood = sum(log(p));
    if abs(log_likelihood-old_log_likelihood)<1e-10
        break
    end
    old_log_likelihood = log_likelihood;
    %% E step
     r = [];
    for i = 1:cluster
        r(:,end+1) = p(:,i)./p_x;
    end
    %% M step
    Nk=[];
    for i = 1:cluster
        Nk(end+1,:) = sum(r(:,i),1);
    end
    u = (sum(r.*pixel))'./Nk;
    Cov = {};
    for i = 1:cluster
        square = (pixel-u(i,:)).*(pixel-u(i,:));
        temp = 0;
        for k = 1:length(pixel)
            temp = temp+r(k,i)*square(k);
        end
        Cov{end+1} = temp/Nk(i,:);
%         Cov{end+1} = r()*((pixel-u(i,:)).*(pixel-u(i,:)))./Nk(i,:);
    end
end
%% 
close all;
seg_pic = zeros(length(idx),1);
for i=1:length(idx)
    if(r(i,1)>r(i,2))
        seg_pic(i) = 0;
    else
        seg_pic(i) = 1;
    end
end
seg_pic = reshape(seg_pic,[h,w,c]);
figure()
subplot(2,2,3)
imshow(seg_pic);
title('GMM segment')
subplot(2,2,1)
imshow(img);
title('original')


subplot(2,2,2)
imshow(kmeansAns/max(max(kmeansAns)))
title('K means segment')

error=kmeansAns/max(max(kmeansAns))-seg_pic/max(max(seg_pic));
subplot(2,2,4)

imshow(error)
title('difference bt k and GMM')

