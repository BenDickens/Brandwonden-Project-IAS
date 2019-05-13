% clear
% nx = 20;
% ny = 40;
% n = nx*ny;
% 
% [Y, X] = meshgrid(linspace(0,2,ny), linspace(-.5,.5,nx));
% x = [X(:)'; Y(:)'];
load x
n = size(x,2);
Omega = zeros(2,2,n);
% cv = [2;1];
% Omega = repmat(cv*cv', [1,1,n]);

for i = 1:n
    b = ([0; 1] - x(:,i))/norm([0; 1] - x(:,i));
    Omega(:,:,i) = b*b';
end

dt = 1e-1;
time = 0:dt:20;
m = size(time,2);

% a = [.25 .75; .6 1];
% Ft = repmat(a'*a, [1,1,n]);
Ft = zeros(2,2,n);
a = rand(2,2,n) - .5;
for i = 1:n
    Ft(:,:,i) = a(:,:,i)'*a(:,:,i);
end

Fd = zeros(n,m);
Fv = zeros(2,n,m);

for i =1:m
    for j = 1:n
        [Dv, Fe] = eig(Ft(:,:,j));
        Fv(:,j,i) = Dv(:, diag(Fe)' == max(max(Fe)));
    end
    itrace = Ft(1,1,:) + Ft(2,2,:);
    Fd(:,i) = itrace(:);
    
    Ft = Ft + 1e-1*Omega*dt;
end
% v = rand(2,n) - .5;
% v(:,max(abs(v))>0) = v(:,max(abs(v))>0) * diag(1./sqrt(sum(v(:,max(abs(v))>0).^2)));

%% Plot
figure(1)
for i = 1:1:m
    cla;
    quiver(x(1,:), x(2,:), Fv(1,:,i), Fv(2,:,i), .5, 'k-', 'ShowArrowHead', 'off');
    hold on
    quiver(x(1,:), x(2,:), -Fv(1,:,i), -Fv(2,:,i), .5, 'k-', 'ShowArrowHead', 'off');
    hold off
    axis([-.5, .5, 0, 2], 'image')
    getframe;
%     pause
end
% interp2(x(1,:), y(2,:), 1);