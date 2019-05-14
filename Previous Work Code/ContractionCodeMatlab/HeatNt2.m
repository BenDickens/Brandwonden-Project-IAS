% du/dt = div(u) + s
% No-flux on northern boundary
% du/dn = u on other boundaries
% Time-dependent source

global p t s np nt delta

xb = [-1.5 1.5];
yb = [-2 2];

%% Define Mesh
load x.mat
p = x';
t = tr;

%% Source function
% s = @(x) (1 - x(1)^2)*(x(2)/2 - x(2)^2/16); %2D-polynomial
% s = @(x) max((1 - 2*norm(x-[0,3]))*3, 0) + max((1 - 2*norm(x-[-.5,1.5]))*3, 0) + max((1 - 2*norm(x-[.5,1.5]))*3, 0); %3 point sources
s = @(x,t) 2*max((1 - 2*norm(x - [0.5*cos(0.5*t), sin(0.5*t)]))*3, 0); %1 point source moving in a circle

%% Triangular elements
nt = size(t,1); %number of triangular elements
np = size(p,1); %number of points
M = zeros(np); %Temporal coefficient matrix
S = zeros(np); %Spatial coefficient matrix
f = zeros(np, 1); %Source vector

delta = zeros(nt,1); %twice the area of each triangle
phi = zeros(2,3,nt); %Test function coefficient matrix
Me = zeros(3,3,nt); %Temporal element matrices
Se = zeros(3,3,nt); %Spatial element matrices
fe = zeros(3,nt); %Element vectors

for i = 1:nt
    %Calculate delta
    delta(i) = (p(t(i,2), 1) - p(t(i,1),1))*(p(t(i,3), 2) - p(t(i,2), 2))-...
        (p(t(i,2), 2) - p(t(i,1), 2))*(p(t(i,3), 1) - p(t(i,2), 1));
    %Create coefficient matrices
    phi(:,:,i) = [circshift(p(t(i,:), 2), -1)' - circshift(p(t(i,:), 2), 1)';...
        circshift(p(t(i,:), 1), 1)' - circshift(p(t(i,:), 1), -1)']/delta(i);
    %Create spatial element matrices
    Se(:,:,i) = -abs(delta(i))/2 * phi(:,:,i)'*phi(:,:,i);
    %Assembly of S
    S(t(i,:), t(i,:)) = S(t(i,:), t(i,:)) + Se(:,:,i); 
    %Create temporal element matrices
    Me(:,:,i) = abs(delta(i))/24 * (ones(3) + eye(3));
    %Assembly of A
    M(t(i,:), t(i,:)) = M(t(i,:), t(i,:)) + Me(:,:,i); 
end

%% Line elements
%Find boundary points
ge = find(abs(p(:,1)-xb(2)) <= 1e-6); gn = find(abs(p(:,2)-yb(2)) <= 1e-6);
gw = find(abs(p(:,1)-xb(1)) <= 1e-6); gs = find(abs(p(:,2)-yb(1)) <= 1e-6);

%Sort according to x or y coordinate
[~, ie] = sort(p(ge, 2)); [~, in] = sort(p(gn, 1));
[~, iw] = sort(p(gw, 2)); [~, is] = sort(p(gs, 1));

%Overwrite saved order
ge = ge(ie); gn = gn(in); gw = gw(iw); gs = gs(is);
l = [ge(1:end-1) ge(2:end);...
    %gn(1:end-1) gn(2:end);...
    gw(1:end-1) gw(2:end);...
    gs(1:end-1) gs(2:end)];

nl = size(l,1);
Sl = zeros(2, 2, np); %Element matrix

for i = 1:nl
    %Create element matrices
    Sl(:,:,i) = -norm(p(l(i,1),:) - p(l(i,2),:))/6 *[2 1; 1 2]; %u' = -u
    %Assembly
    S(l(i,:), l(i,:)) = S(l(i,:), l(i,:)) + Sl(:,:,i); 
end

%% Temporal integration
dt = 1e-1; %time step
T = 300; %Number of time steps
ti = 0:dt:(T-1)*dt; %Time vector

u = zeros(np, T); %Solution

f_old = AssembleVector(ti(1));

for i = 1:T-1
    f_new = AssembleVector(ti(i+1));
%     u(:,i+1) = M\( (M + dt*S)*u(:,i) + dt*f_old); %Euler Forward
%     u(:,i+1) = (M - dt*S)\(M*u(:,i) + dt*f_new); %Euler Backward
    u(:,i+1) = (M - .5*dt*S)\((M + .5*dt*S)*u(:,i) + dt*0.5*(f_old + f_new)); %Crank-Nicholson
    f_old = f_new;
end

%% Plot source function
% figure(2)
% set(figure(2), 'name', 'Source Function')
% sp = zeros(np,1);
% for i = 1:np
%     sp(i) = s(p(i,:));
% end
% trisurf(t, p(:,1), p(:,2), sp)
% axis equal

%% Plot solution
figure(1)
set(figure(1), 'name', 'Solution')
for i = 1:T
    clf    
    axis([-1 1 0 4 0 1], 'image')
    trisurf(t, p(:,1), p(:,2), u(:,i))
    view(3)
    getframe;
end