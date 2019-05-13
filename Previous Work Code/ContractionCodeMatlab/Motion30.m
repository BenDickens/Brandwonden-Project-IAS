clear
tic
%% Input
r = 1e-2; %Radius
n = 100; %Number of cells
% cx(:,:,1) = [0.5*(-1).^(1:n/2), 0.5*linspace(-1+2*r,1-2*r,n/2); linspace(2*r,2-2*r,n/2), 1+(-1).^(1:n/2)]; %Initial distribution
cx(:,:,1) = [0.5*(-1).^(1:n/2), linspace(-.5+2*r,.5-2*r,n/2); linspace(2*r,2-2*r,n/2), 1+(-1).^(1:n/2)]; %Initial distribution
dt = 1e-3; %Time step

%% Movement
for i = 1:3000
%     cv = [0.1*sin(i/1000 * 2 * pi); 1+0.1*cos(i/1000 * 2 * pi)]*ones(1,n) - cx(:,:,i); %Moving origin
    cv = [0;1]*ones(1,n) - cx(:,:,i); %Direction to origin
    
    cv(:, max(abs(cv)) > 0) = cv(:, max(abs(cv)) > 0) * diag(1./sqrt(sum(cv(:, max(abs(cv)) > 0).^2))); %Normalize cv
    cv = cv + 1*(randn(2,n)); %Add random deviation
    cv(:, max(abs(cv)) > 0) = cv(:, max(abs(cv)) > 0) * diag(1./sqrt(sum(cv(:, max(abs(cv)) > 0).^2))); %Normalize cv
    
    D = dist(cx(:,:,i)); %Distance matrix
    [c.de, c.of] = find((D <= 2*r + 2*dt).*(D>0)); %Find colliding cells
    c.in = find((D <= 2*r + 2*dt).*(D>0)); %Linear indexing
    
    if ~isempty(c.in)
        di = zeros(2,n^2);
        di(:,c.in) = cx(:,c.de,i) - cx(:, c.of,i); %get vector between cells
        
        theta = zeros(n);
        theta(c.in) = atan2(di(2,c.in), di(1,c.in)) - atan2(cv(2,c.of), cv(1,c.of)); %Calculate theta for each offending cell
        theta(abs(theta) > pi) = theta(abs(theta) > pi) - sign(theta(abs(theta) > pi))*2*pi; %Modulo 2pi
        theta(abs(theta)> 0.5*pi) = NaN; %An obtuse collision is not a collision
        
        theta(:, c.of(abs(theta(c.in)) == 0)) = 0; %Head-on collision makes theta zero
        
        [~, I] = max(abs(theta));
        thetam = theta(I + (0:n:n^2-1))'; %Save maximum theta so minimum deflection can be chosen
        di = di(:,I + (0:n:n^2-1)); %Find corresponding cell-cell vector
        
        thetam(max(theta).*min(theta) < 0) = 0; %Blocked by two cells
        
        c.of = setdiff(c.of, find(thetam'==0 & sum(isnan(theta)) >=1)); %Take out obtuse collisions
        thetam = thetam(c.of); %Only offending cells are interesting
        di = di(:,c.of); %Idem
        
        A = zeros(2*length(c.of), 2); %Create rotation vector
        
        A(2:2:end,1) = -sign(thetam);
        A(1:2:end,2) = sign(thetam);
        C = A*di; %Rotate cell-cell vectors according to thetam
                      
        if length(c.of)==1;
            cv(:,c.of) = C(:,1); %Bypass diag compatibility problem
        elseif length(c.of)>1;
            cv(:, c.of) = [diag(C(1:2:end,:), 0)'; diag(C(2:2:end,:), 0)']; %Create new velocity vectors
        end
        
        cv(:, max(abs(cv)) > 0) = cv(:, max(abs(cv)) > 0) * diag(1./sqrt(sum(cv(:, max(abs(cv)) > 0).^2))); %Normalize cv
        
        D = dist(cx(:,:,i) + dt*cv); %Final check for collisions in next step
        [c.de, c.of] = find((D < 2*r).*(D>0));
        cv(:, c.of) = 0; %Stop colliding cells
        
    end
    
    cx(:,:,i+1) = cx(:,:,i) + dt*cv; %Create new coordinates
end
toc
%% Plot
figure(1)
for i = 1:10:3000
    scatter(cx(1,:,i), cx(2,:,i));
    axis([-2, 2, -2, 2]);
    getframe;
end