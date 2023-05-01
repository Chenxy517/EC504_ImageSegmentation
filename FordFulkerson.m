% function [f, cap] = FordFulkerson(CapMatrix, S, T)

% CapMatrix = [0,2,0,100;0,0,0,100;50 100 0 0;0 0 0 0];
% S = 3;T=4;
% CapMatrix = ...
% [0 70 5 0 100 1;...
% 70 0 0 8 90 2;...
% 5 0 0 80 1 90;...
% 0 8 80 0 2 80;...
% 100 90 1 2 100 0;...
% 1 2 90 80 0 100];
% S = 5;T=6;
[h,w,c] = size(img);
CapMatrix = graph;
graph = CapMatrix;
S = length(graph)-1;
T = length(graph);
cap = CapMatrix;
f = 0;

while true
    
    path = findPath(cap, S, T);
    

    if path(1) == 0
        break
    end
    
    flow = max(max(cap));
    
    for j = 2:length(path)
        flow = min(flow, cap(path(j), path(j-1)));
    end
    
    for j = 2:length(path)
        a = path(j); b = path(j-1);
        cap(a,b) = cap(a,b) - flow;
        cap(b,a) = cap(b,a) + flow;
    end
    
    f = f + flow;

end
cut = [];
% 
for i=1:length(cap)-2
    for j =1:length(cap)-2
        if graph(i,j)~=0&& cap(i,j)==0
            cut = [cut;i;j];
        end
    end
end

boundary = zeros(length(cut),2);
for i = 1:length(cut)
    boundary(i,:) = [floor(cut(i)/w)+1,cut(i)-floor(cut(i)/w)*w];
end
[hh,ww] = size(boundary);
boundary_inv = boundary;%zeros(size(boundary));
for i = 1:hh
    for j = 1:ww
        boundary_inv(i,2) = 60-boundary(i,2);
    end
end
figure()
imshow(img)
hold on
scatter(boundary(:,1),boundary(:,2),5,'filled');
% end

function F = findPath(CapMatrix, S, T)
% BFS (Breadth-first Search)

n = length(CapMatrix);

front = 1; 
back = 2;

queue = zeros(1,n);
queue(front) = S;

pred = zeros(1,n);
pred(S) = S; 

while front ~= back
    node = queue(front);
	front = front + 1;
    
    for i = 1:n
        if pred(i) == 0 && CapMatrix(node,i) > 0
            queue(back) = i;
            back = back + 1;
            pred(i) = node;
        end
    end


end
path = zeros(1,n);

if pred(T) ~= 0
    i = T; c = 1;
	
    while pred(i) ~= i
    	path(c) = i;
        c = c + 1;
        i = pred(i);
    end
    
    path(c) = S;
    path(c+1:n) = [];
end
F = path;
end
