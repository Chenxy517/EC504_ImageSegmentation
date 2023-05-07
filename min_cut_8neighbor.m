clc; clear all;close all;
lambda = 0.5;
SIGMA = 0.1;
cluster = 2;

% img = im2double(imread("C:\Users\Yang Shu\Desktop\tree.jpg"));
img = load('img.mat');
img = struct2array(img);
img = imresize(img, [50 50]);
img = img(:,:,1);
[h,w,c] = size(img);
pixel = reshape(img',[h*w,c]);
mean = mean(pixel);
std = std(pixel);
pixel = (pixel-mean)./std;
pixel = (pixel-mean);
[idx, class] = kmeans(pixel, cluster);
kmeansAns = reshape(idx,[h,w]);
[u,Cov,r] = GMM_seg(idx, class,pixel,img,lambda);
pixel_borad = [pixel ; u];
graph = build_graph(pixel_borad,h,w,c,SIGMA,r,lambda);

function graph = build_graph(pixel_borad,h,w,c,SIGMA,r,lambda)
    graph = zeros(length(pixel_borad));
    graph = graph;
    % calculate the N link cost
    K = -inf;
    for i = 1:h
        for j =1:w

            if i<h
                ec = edge_cost(pixel_borad((i-1)*w+j),pixel_borad(i*w+j),SIGMA);
                graph((i-1)*w+j,i*w+j) = ec;graph(i*w+j,(i-1)*w+j)=ec;
                K = max(K, ec);
            end
            if j<w
                ec = edge_cost(pixel_borad((i-1)*w+j),pixel_borad((i-1)*w+j+1),SIGMA);
                graph((i-1)*w+j,(i-1)*w+j+1) = ec;graph((i-1)*w+j+1,(i-1)*w+j) = ec;
                K = max(K, ec);
            end
            if i<h && j<w
                ec = edge_cost(pixel_borad((i-1)*w+j),pixel_borad(i*w+j+1),SIGMA);
                graph((i-1)*w+j,i*w+j+1) = ec/sqrt(2);graph(i*w+j+1,(i-1)*w+j) = ec/sqrt(2);
                K = max(K, ec);
            end
            graph((i-1)*w+j,(i-1)*w+j) = 0;
        end
    end
    % calculate the T link cost
    r_lambda = r*lambda;
    for i = 1:length(pixel_borad)-2
            graph(i,end-1) = r_lambda(i,1);
            graph(i,end) = r_lambda(i,2);
            graph(end-1,i) = r_lambda(i,1);
            graph(end,i) = r_lambda(i,2);
    end
    graph(end-1,end-1)=K;
    graph(end-1,end)=0;
    graph(end,end-1)=0;
    graph(end,end)=K;
end

function ec = edge_cost(ip,iq,SIGMA)
    ec = 3 * exp(- (ip - iq )^2 / (2 * SIGMA^2));
end

% function x = GMM_seg(cluster,idx, class,pixel)
%     int_cov,initial_priors = init_GMM(cluster,pixel);
%     
%     x=1;
% end
% 
% function int_cov,initial_priors = init_GMM(cluster,pixel)
%     int_cov = {};
%     initial_priors = zeros(cluster,1);
%     for i = 1:cluster
%         data = [];
%         for j=1:size(idx)
%             if idx(j) == i
%                 data(end+1,:)=pixel(j,:);
%             end
%         end
%         int_cov{end+1} = cov(data);%+lambda*diag(ones(1, c));
%         initial_priors(i) = length(data)/length(idx);
%     end
% end

function [u,Cov,r] = GMM_seg(idx, class,pixel,img,lambda)
%     lambda = 1*10^-5;
%     img = im2double(imread("C:\Users\Yang Shu\Desktop\cat1.png"));
%     img = imresize(img, [50 50]);
    cluster = 2;
    %figure()
    % imshow(img);
    img = img(:,:,1);
    [h,w,c] = size(img);
    pixel = reshape(img,[h*w,c]);
    Mean = mean(pixel);
    Std = std(pixel);
    pixel = (pixel-Mean)./Std;
    pixel = (pixel-Mean);
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
    seg_pic = zeros(length(idx),1);
    for i=1:length(idx)
        if(r(i,1)>r(i,2))
            seg_pic(i) = 0;
        else
            seg_pic(i) = 1;
        end
    end
    seg_pic = reshape(seg_pic,[h,w,c]);
end