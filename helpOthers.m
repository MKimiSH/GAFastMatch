RGB = imread('1.jpg');

[y x] = find(RGB > 100);


% x0 = min(x)
% x1 = max(x)
% y0 = min(y)
% y1 = max(y)
% 
% x = x - mean(x);
% x = x / std(x);
% y = y - mean(y);
% y = y / std(y);

x0 = min(x)
x1 = max(x)
y0 = min(y)
y1 = max(y)

X = zeros(size(x,1),2);
X(:,1) = x;
X(:,2) = y;

F=@(p,X)(p(1)*X(:,1).^2+p(2)*X(:,2).^2+p(3)*X(:,1)+p(4)*X(:,2)+p(5));
p0 = [1 1 1 1 10]
p=nlinfit(X,zeros(size(X,1),1),F,p0);

figure;
h = ezplot(@(X,y)F(p,[X,y]),[-1+x0,x1+1,-1+y0,y1+1]);

set(h,'color','y','linewidth',3);