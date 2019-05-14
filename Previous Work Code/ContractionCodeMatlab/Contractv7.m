%% Model with contraction
% Added visco-elastic relationship

%% Constants
%% Time and Space
dt = 5e-1; %time step
T = floor(10*24/dt) + 1; %Number of time steps
ti = 0:dt:(T-1)*dt; %Time vector

xb = [-3 3];
yb = [-2 2];

%% Cumming
DPDGF = 2000/(400^2); %Diffusion of PDGF

DtPA.max = 4000/(400^2);  %Maximum diffusion tPA 1000
DtPA.min = 200/(400^2);   %Minimum diffusion tPA 50
DTGF.max = 3000/(400^2);  %Maximum diffusion TGF 750
DTGF.min = 300/(400^2);   %Minimum diffusion TGF 75

vmmax = 400/400; %Maximum velocity macrophages
vfmax = 150/400; %Maximum velocity fibroblasts
vfi = 0.4; %Factor for influence of fibrin in cell velocity 0.4

fg = 2; %Directional sensitivity to gradient fibroblasts
sqrtalpha = sqrt(1e-2); %Deviation

rfc1 = 400/(400^2); %Minimum collagen synthesis (30/(400^2))
rfc2 = 800/(400^2); %Maximum collagen synthesis (60/(400^2))

spa = 5; %Probability factor for macrophage spawn

rtgf = 150/(400^2); %Maximum excretion of TGF 5/(400^2)

sf = 20/400*0.1; %Source function tPA

mit_n = .5; %Mitosis receptor threshold 
mit_t = 48; %Mitosis time threshold

rp = 0.15; %fibrin decay due to tPA

n_b = .1; %TGF expression due to fibrin concentration .1
n_bc = 1; %Rate of unbinding receptors
n_gc = 5; %Rate of binding receptors

%% Contraction constants
umax = 0; %Maximum stress over all time steps
E = 16e-3; %Elasticity modulus 100kPa
nu = .3; %Poisson's ratio

mu = E/(2 + 2*nu); % Shear modulus
K =  E/(3*(1 - 2*nu)); % Bulk modulus
mu_1 = 2*mu;
mu_2 = K - 2*mu/3;

fc.min = 1e-6 * 40^2 / pi; % Minimum force per fibroblast
fc.max = 3e-6 * 40^2 / pi; % Maximum force per myofibroblast

srho = 1; %Spring constant boundary condition

tau_max = 12; % Maximum time for plastic forces
fpl_2 = 10; % Factor plastic forces / temp forces

%% Create initial distribution of cells
rc = 10/400; %Radius cells

nf = 250; %Total possible number of fibroblasts
nfi = 75; %Initial number of fibroblasts
nm = 100; %Total possible number of macrophages
nc = nm + nf; %Total number of cells

c.do = 1; %Cell indices outside domain
c.co = 1; %Colliding cell indices

xcf = zeros(nfi, 2, T); %Preallocate cell positions

%Create initial distribution
xcf(:,1,1) = 6*(rand(nfi, 1)-.5);
xcf(:,2,1) = 4*(rand(nfi, 1)-.5); %Homogeneous
% xc(nm+1:nm+nfi,2,1) = (2 -4*sqrt(1 - rand(nfi, 1))); %Linear

fprintf('Distributing cells... ');
while ~isempty(c.co) || ~isempty(c.do)
    Di = dist(xcf(:,:,1)'); %Distance matrix
    [c.co, ~] = find((Di <= 2*rc).*(Di>0)); %Find colliding cells
    [c.do, ~] = find((abs(xcf(:,1,1))<=2).*(xcf(:,2,1) >= 0)); %Find cells in wound space
    
    xcf([c.do; c.co],1,1) = 6*(rand(length(c.do)+length(c.co),1)-.5);
    xcf([c.do; c.co],2,1) = 4*(rand(length(c.do)+length(c.co),1)-.5);
end

xc = [NaN(nm, 2, T); xcf; NaN(nf - nfi, 2, T)];

% Vector for myo-status
myo = zeros(nf, 1);
myo(1:nfi) = .5 - .25*xcf(:,2); % Linear increase in depth

fprintf('done\n');

%% Define Mesh
fprintf('Creating mesh... ');

try
    isempty(tr); %Only if mesh does not exist yet
    
catch err
    yb = yb + 4; %Coordinate transform for mesh creation
    xh = 0.1; %Minimum distance in mesh
    fd = @(p) drectangle0(p,xb(1),xb(2),yb(1),yb(2));
    fh = @(p) min(.5*xh + .4*xh*abs(dpoly(p,[-2, yb(2); -2, mean(yb); 2, mean(yb); 2, yb(2); -2, yb(2)])), xh);
    %         xh+xh*abs(dpoly(p,[xb(1), yb(2); xb(1), yb(1); xb(2), yb(1); xb(2), yb(2)]))), xh); %Refinement on boundaries
    figure(4)
    [x,tr] = distmesh2d(fd, fh, xh, [xb; yb],[]);
    
    x(:,2) = x(:,2) - 4;
    yb = yb - 4;
    
    wound = [-2, yb(2); -2, mean(yb); 2, mean(yb); 2, yb(2)];
    hold on
    plot(wound(:,1), wound(:,2) + 4, 'r');
    hold off
end

ntr = size(tr, 1);
np = size(x, 1);

xp = zeros(np, 2, T); %Moving mesh coordinates
xp(:,:,1) = x;

fprintf('done\n');
tic

%% Triangular element index lists
Imat = [1 1 1;2 2 2;3 3 3];
Jmat = Imat';
Itopo = tr(:,Imat);
Jtopo = tr(:,Jmat);

%% Line elements
%Find boundary points
ge = find(abs(x(:,1)-xb(2)) <= 1e-6); gn = find(abs(x(:,2)-yb(2)) <= 1e-6);
gw = find(abs(x(:,1)-xb(1)) <= 1e-6); gs = find(abs(x(:,2)-yb(1)) <= 1e-6);

%Sort according to x or y coordinate
[~, ie] = sort(x(ge, 2)); [~, in] = sort(x(gn, 1));
[~, iw] = sort(x(gw, 2)); [~, is] = sort(x(gs, 1));

%Overwrite saved order
ge = ge(ie); gn = gn(in);
gw = gw(iw); gs = gs(is);
l = [ge(1:end-1), ge(2:end);...
    gw(1:end-1), gw(2:end);...
    gs(1:end-1), gs(2:end)];

ln = [gn(1:end-1), gn(2:end)];

trn = zeros(length(gn) - 1, 1); %Find northern boundary triangular elements
for e = 1:length(gn) - 1
    trn(e) = find(any(tr' == gn(e))' .* any(tr' == gn(e+1))');
end

nl = size(l,1); %Number of line elements

Iltopo = l(:, [1 1;2 2]);
Jltopo = l(:, [1 1;2 2]');

%% Preallocation
PDGF = zeros(np, T); %Solution PDGF
tPA = zeros(np, T); %Solution tPA
TGF = zeros(np, T); %Solution TGF
fi = zeros(np, T); %Solution Fibrin
collas = zeros(np, T); %Solution Collagen density
collat = zeros(np, 4, T); %Solution Collagen tensor

nmPDGF = zeros(nm, T); %Solution receptor bindings macrophages
nftgf = zeros(nf, T); %Solution receptor bindings fibroblasts
mit_f = zeros(nf, 1); %Amount of time bounded receptors are above threshold

StPA = zeros(np, 1); %Integral over tPA from 0 to t

vm = zeros(nm,2); %Velocity macrophages
vf = zeros(nf,2); %Velocity fibroblasts

%Element matrices
Mlist = zeros(ntr,9); %Preallocate element matrix for M
phi = zeros(3,3,ntr); %Test function coefficients matrix

Flist = zeros(ntr, 9); %tPA
Fl = zeros(nl, 4); %Line elements tPA
Sclist = zeros(ntr,9);
Llist = zeros(nl, 4);

xxlist = zeros(ntr, 9);
yylist = zeros(ntr, 9);
xylist = zeros(ntr, 9);
Gslist = zeros(nl, 4);

reverseStr = ''; %For status update

%% Create mass and coefficient matrices
X = xp(:,1,1);
Y = xp(:,2,1);

% Build matrices
Xtri = X(tr);
Ytri = Y(tr);
delta = 2*polyarea(Xtri,Ytri,2);

for e = 1:ntr
    %Create coefficient matrices
    Q = [ones(3,1), Xtri(e,:)', Ytri(e,:)'];
    phi(:,:,e) = Q\eye(3);
    %Create temporal element matrices
    %     Me = delta(e)/24 * (ones(3) + eye(3)); %Holand Bell
    Me = delta(e)/6 * eye(3); %Newton-Cotes
    %Assembly of M
    Mlist(e,:) = Me(:);
end

M_old = sparse(Itopo, Jtopo, Mlist, np, np); %Mass matrix

%% Initial conditions
PDGF(:,1) = .5*(1 + tanh((2-X)/.2))*.5.*(1 + tanh((2+X)/.2))*.5.*(1 + tanh(Y/.2));
fi(:,1)   = .5*(1 + tanh((2-X)/.2))*.5.*(1 + tanh((2+X)/.2))*.5.*(1 + tanh(Y/.2));
collas(:,1) = 1 - fi(:,1);
collat(:, [1, 4], 1) = [collas(:,1)/2, collas(:,1)/2];
u_old = zeros(np,2);
g_pl = zeros(np,2);
tau = zeros(np,1);

%% Define wound area
nw = 150; 
xw = zeros(nw,2,T);
wt = linspace(0,1,nw)';
xw(:,:,1) = [(wt<=.25)*wound(1,1) + (abs(wt - .5) < .25).*(8*wt - 4) + (wt>=.75)*wound(1,2), ...
    (wt<=.25).*(2 - 8*wt) + (wt>=.75).*(8*wt - 6)];

xwe = zeros(nw, 1);
for e = 1:ntr %Indexing wound points to corresponding elements
    eiw = find(inpolygon(xw(:,1,1), xw(:,2,1), Xtri(e,:), Ytri(e,:))); %Find which points are in the element
    xwe(eiw) = e;
end

phi_w = zeros(nw, 3);
for e = 1:nw
    phi_w(e,:) = [1, xw(e ,:,1)]*phi(:,:,xwe(e));
end

X2w = sparse( [1:nw; 1:nw; 1:nw]', tr(xwe,:), phi_w, nw, np);

%% Temporal Integration
fprintf('Temporal integration started at '); 
disp(datestr(now));

for i = 1:T-1
    %% Matrices for interpolation and contraction
    
    xce = zeros(nc, 1); %List of elements for each cell
    
    for e = 1:ntr
        eic = find(inpolygon(xc(:,1,i), xc(:,2,i), Xtri(e,:), Ytri(e,:))); %find which cells are in the element
        xce(eic) = e;
        
        phixx = delta(e)/2 * phi(2,:,e)'*phi(2,:,e);
        xxlist(e,:) = phixx(:);
        phiyy = delta(e)/2 * phi(3,:,e)'*phi(3,:,e);
        yylist(e,:) = phiyy(:);
        phixy = delta(e)/2 * phi(2,:,e)'*phi(3,:,e);
        xylist(e,:) = phixy(:);
    end

    Diffxx = sparse(Itopo, Jtopo, xxlist ,np, np);
    Diffyy = sparse(Itopo, Jtopo, yylist ,np, np);
    Diffxy = sparse(Itopo, Jtopo, xylist ,np, np);
    Diffyx = Diffxy';
    
    c_xce = find(xce ~=0); %Existing cell indices
    
    phi_c = zeros(length(c_xce),3); %Preallocation
    gradphi1_c = zeros(length(c_xce),3);
    gradphi2_c = zeros(length(c_xce),3);
    flist = zeros(length(c_xce(c_xce > nm)),6);
    gpllist = zeros(length(c_xce(c_xce > nm)),6);
    
    for e = 1:length(c_xce)
        phi_c(e,:) = [1, xc(c_xce(e) ,:, i)]*phi(:,:,xce(c_xce(e)));
        gradphi1_c(e,:) = phi(2,:,xce(c_xce(e)));
        gradphi2_c(e,:) = phi(3,:,xce(c_xce(e)));
    end
    
    X2C = sparse([c_xce, c_xce, c_xce], tr(xce(c_xce), :), phi_c, nc, np);
    Grad1 = sparse([c_xce, c_xce, c_xce], tr(xce(c_xce), :), gradphi1_c, nc, np);
    Grad2 = sparse([c_xce, c_xce, c_xce], tr(xce(c_xce), :), gradphi2_c, nc, np);
    
    %% Interpolation
    
    fi_c = X2C*fi(:, i); % Fibrin values
    colla_c = X2C*collas(:, i); % Collagen density values
    collat_c = X2C*collat(:,:, i); % Collagen tensor        
    PDGF_c = X2C*PDGF(:, i); % PDGF values
    
    gradPDGF_c = [Grad1(1:nm, :)*PDGF(:, i), Grad2(1:nm, :)*PDGF(:,i)]; %Gradient of PDGF
    gradTGF_c = [Grad1(nm+1:end, :)*TGF(:, i), Grad2(nm+1:end, :)*TGF(:,i)]; %Gradient of TGF
    
    gradPDGF = gradPDGF_c./((1 + sqrt(gradPDGF_c(:,1).^2 + gradPDGF_c(:,2).^2))*ones(1,2));
    gradTGF = gradTGF_c./((1 + sqrt(gradTGF_c(:,1).^2 + gradTGF_c(:,2).^2))*ones(1,2));
    
    %% Calculate contraction
    f_xce = c_xce(c_xce > nm); %Find existing fibroblasts
    
    %Calculate strain on cell positions
    ux = Grad1*u_old(:,1);
    vy = Grad2*u_old(:,2);
    
    fce = myo*fc.max + (1 - myo)*fc.min;
    
    for e = 1:length(f_xce)
        dxc = ones(3,1)*xc(f_xce(e),:, i) - xp(tr(xce(f_xce(e)), :), :, i); %Find direction of force
        r = sqrt(dxc(:,1).^2 + dxc(:,2).^2); %Distance to mesh points
        
        phi_f = [1, xc(f_xce(e) ,:, i)]*phi(:,:,xce(f_xce(e)));
        f1 = max(1,colla_c(f_xce(e))) * fce(f_xce(e)) * diag(phi_f)*diag(1./r)*dxc;
        flist(e, :) = f1(:);
        gpllist(e, :) = (ux(f_xce(e)) + vy(f_xce(e)) >0) * f1(:); %Knots are only formed if there is a stress on the tissue
    end
    
    umax = max(umax, max(ux + vy)); %Save maximum stress for analysis
    g = sparse([tr(xce(f_xce), :), tr(xce(f_xce), :)], ones(length(f_xce), 1)*[1, 1, 1, 2, 2, 2], flist, np, 2); %Contractile forces
    gpl = sparse([tr(xce(f_xce), :), tr(xce(f_xce), :)], ones(length(f_xce), 1)*[1, 1, 1, 2, 2, 2], gpllist, np, 2); %Plastic forces
        
    g_pl = g_pl + fpl_2*dt/tau_max * diag(tau<tau_max) * gpl;
    tau = tau + (sum(abs(gpl),2) > 0)*dt;
    
    for e = 1:nl
        %Line elements for u
        Gsl = srho * norm(xp(l(e,1),:,i) - xp(l(e,2),:,i)) / 6 * (ones(2)+eye(2));
        Gslist(e,:) = Gsl(:);
    end
    
    Gs = sparse(Iltopo, Jltopo, Gslist, np, np);
    
    %Spatial matrices
    S1 = E/((1 + nu)*(1 - 2*nu))*((1 - nu)*Diffxx + .5*(1 - 2*nu)*Diffyy) + Gs;
    S2 = E/((1 + nu)*(1 - 2*nu))*(nu*Diffxy + .5*(1 - 2*nu)*Diffyx);
    S3 = E/((1 + nu)*(1 - 2*nu))*(nu*Diffyx + .5*(1 - 2*nu)*Diffxy);
    S4 = E/((1 + nu)*(1 - 2*nu))*((1 - nu)*Diffyy + .5*(1 - 2*nu)*Diffxx) + Gs;
    S = [S1, S2; S3, S4];
    
    %Temporal matrices
    T1 = (mu_1 + mu_2)*Diffxx + .5*mu_1*Diffyy;
    T2 = mu_2*Diffxy + .5*mu_1*Diffyx;
    T3 = mu_2*Diffyx + .5*mu_1*Diffxy;
    T4 = (mu_1 + mu_2)*Diffyy + .5*mu_1*Diffxx;
    TM = [T1, T2; T3, T4];
    
    %Solve for u and reshape
    u2 = (TM + dt*S)\(TM*u_old(:) + dt*(g(:) + g_pl(:)));
    u = reshape(u2, np, 2);
    
    u_delta = u - u_old; %Change in contraction compared to last iterate
    
    %Update mesh
    xp(:,:,i+1) = xp(:,:,1) + u;
    
    %% Cell movement
    vms = vmmax *(5/8 - 3/16*xc(1:nm,2,i)).*(1 - vfi * fi_c(1:nm)).* (.25 + 3./(4 + 80*(1-nmPDGF(:,i)).^6)); %Velocity size
    vfs = vfmax *(5/8 - 3/16*xc(nm+1:end,2,i)).*(1 - vfi * fi_c(nm+1:end)).* (.25 + 3./(4 + 80*(1-nftgf(:,i)).^6)); %Velocity size
    
    vm = diag(1 - nmPDGF(:,i))*gradPDGF + sqrtalpha * randn(nm,2); %Direction vector
    vf = fg*diag(1 - nftgf(:,i))*gradTGF  +  sqrtalpha * randn(nf,2);
    
    vfc = zeros(nf, 2); %Velocity vector after contact guidance
    for e = 1:nf
        %Add contact guidance for fibroblasts
        Omegahat = reshape(collat_c(nm + e,:), 2,2);
        Omegahat2 = Omegahat * diag(1./sqrt(sum(Omegahat.^2)));
        vfc(e,:) = (((1 - colla_c(nm + e))*eye(2) + colla_c(nm + e)*Omegahat2)*vf(e,:)')';
    end
    
    vcs = ([vms; vfs]*ones(1,2)).*[vm ; vfc];
    
    xcp = xc(:,:,i) + dt*vcs + X2C*u_delta; %Predicted new values
    Di = dist(xcp'); %Distance matrix in prediction
    
    SortedDi = sort(Di, 2);
    minDi = SortedDi(nm + 1:end, 2);
    
    [c.de, c.of] = find((Di < 2*rc).*(Di>0)); %Find colliding cells
    c.in = find((Di < 2*rc).*(Di>0)); %Linear indexing
    
    if ~isempty(c.in)
        di = zeros(nc^2,2);
        di(c.in,:) = xcp(c.de,:) - xcp(c.of,:); %get vector between cells        
        def = zeros(length(c.in), 2);
        
        for j = 1:length(c.in)
            def(j,:) = 1/dt*(Di(c.in(j))/(2*rc) - 1)*di(c.in(j),:);
        end
        
        vcs = vcs + sparse([c.of, c.of], ones(length(c.of),1)*[1 2], def, nc, 2);
    end
    
    xcp = xc(:,:,i) + dt*vcs + X2C*u_delta; %Predicted new values
    
    xcn = zeros(nc,1); %Find cells in northern region    
    for e = trn'
        xcn = xcn + (xce == e);
    end
    
    for e = find(xcn)'
        lne = ln(trn == xce(e),:); %Line element corresponding to cell
        lnev = xp(lne(2), :, i+1) - xp(lne(1), :, i+1);
        k = lnev(1) .* (xcp(e,2) - xp(lne(1),2,i+1) + rc) < ...
            lnev(2) .* (xcp(e,1) - xp(lne(1),1,i+1));
        vcs(e,:) = k*vcs(e,:) + ...
            (1-k)*((vcs(e,:)*lnev')*lnev + ...
            rc*[lnev(2), -lnev(1)])/norm(lnev) ; %Stop cells leaving domain
    end

    xc(:,:,i+1) = xc(:,:,i) + dt*vcs + X2C*u_delta; %Calculate new cell coordinates
    
    %% Fibroblast mitosis and death
    mit_f = mit_f + dt *(nftgf(:,i) > mit_n);
    
    c.mit = find((mit_f >= mit_t).*(minDi > 2*rc));
    fnew = length(c.mit);    
    c.nan = find(isnan(xc(nm+1:end, 1, i)));
    
    xc(nm + c.nan(1:fnew), :, i+1) = xc(nm + c.mit, :, i+1);
    myo(nm + c.nan(1:fnew)) = myo(c.mit);
    
    mit_f(c.mit) = 0;
    
    live_f = find(~isnan(xc(nm+1:end, 1, i)));
    death_f = false(nc,1);
    death_f(nm + live_f) = (ux(live_f + nm) + vy(live_f + nm) < -0.5);
    xc(death_f, :, i+1) = NaN;
    
    %% Macrophage spawn and decay
    
    PDGF_cw = X2w*PDGF(:,i);
    
    m_spawn = (binornd(1, spa/nw*dt*PDGF_cw) == 1);
    c.nan = find(isnan(xc(1:nm, 1, i)));
    c.live = find(~isnan(xc(1:nm, 1, i)));
    
    xc(c.nan(1:sum(m_spawn)), :, i+1) = xw(m_spawn, :,i);
    
    death_m = false(nm,1);
    death_m(c.live) = (exprnd(72 * PDGF_c(c.live)) < dt);
    xc(death_m, :, i+1) = NaN;
    
    %% Create new phi and mass matrices
    X = xp(:,1,i+1);
    Y = xp(:,2,i+1);
    
    % build matrices
    Xtri = X(tr);
    Ytri = Y(tr);
    delta = 2*polyarea(Xtri,Ytri,2);
    
    for e = 1:ntr
        %Create coefficient matrices
        Q = [ones(3,1), Xtri(e,:)', Ytri(e,:)'];
        phi(:,:,e) = Q\eye(3);
        %Create temporal element matrices
        %     Me = delta(e)/24 * (ones(3) + eye(3)); %Holand Bell
        Me = delta(e)/6 * eye(3); %Newton-Cotes
        %Assembly of M
        Mlist(e,:) = Me(:);
        %Elements for Sc
        Sce = delta(e)/2 * phi(2:3,:,e)'*phi(2:3,:,e);
        Sclist(e,:) = Sce(:)';
    end
    
    M = sparse(Itopo, Jtopo, Mlist, np, np); %Mass matrix
    
    for e = 1:nl
        %Line elements for Sc
        Scl = norm(xp(l(e,1),:,i+1) - xp(l(e,2),:,i+1)) / 6 * (ones(2)+eye(2)); %u' = -u;
        Llist(e,:) = Scl(:)';
    end
    
    Sc = sparse([Itopo(:); Iltopo(:)], [Jtopo(:); Jltopo(:)], [Sclist(:); Llist(:)], np, np); %Stiffness matrix

    %% Move wound boundary
    phi_w = zeros(nw, 3);
    for e = 1:nw
        phi_w(e,:) = [1, xw(e ,:, i)]*phi(:,:,xwe(e));
    end
    
    X2w = sparse( [1:nw; 1:nw; 1:nw]', tr(xwe,:), phi_w, nw, np);
    xw(:,:,i+1) = xw(:,:,1) + X2w*u;
    
    %% PDGF
    PDGF(:,i+1) = (M + dt*DPDGF*Sc)\(M_old*PDGF(:,i));
    
    %% Picard iteration for tPA and fibrin
    stPA = X2w'*sf*8/nw*ones(nw,1); %Source function tPA (length wound boundary = 8)
    
    %Initial guess is current value
    p.fi = fi(:,i);
    
    stop = 0;
    while ~stop
        for e = 1:ntr
            %Assembly of F
            Fe = delta(e)/6 * phi(2:3,:,e)'*phi(2:3,:,e) * sum(p.fi(tr(e,:)));
            Flist(e,:) = Fe(:)';
        end
        for e = 1:nl
            Fl = norm(x(l(e,1),:) - x(l(e,2),:))/3 * (p.fi(l(e,1))*[3 1; 1 1] + p.fi(l(e,2))*[1 1; 1 3]);
            Llist(e,:) = Fl(:)';
        end
        
        F = sparse([Itopo(:); Iltopo(:)], [Jtopo(:); Jltopo(:)], [Flist(:); Llist(:)], np, np); %Spatial coefficient matrix tPA
        
        %Create next iterates using Euler Backward
        tPA(:,i+1) = (M + dt*DtPA.max*Sc + dt*(DtPA.min - DtPA.max)*F)\(M_old*tPA(:,i) + dt*stPA );
        fi(:,i+1) = fi(:,1).* exp( - rp*dt*(StPA + .5*(tPA(:,i) + tPA(:,i+1)) )); %Semi-analytic solution
        
        stop = (norm(fi(:,i+1) - p.fi)/norm(p.fi) <= 1e-5);
        
        %Save next iterate
        p.fi = fi(:,i+1);
    end
    %Update the integral over tPA
    StPA = StPA + .5*(tPA(:,i) + tPA(:,i+1));
    
    %% New interpolation for cells
    
    xce = zeros(nc, 1);
    for e = 1:ntr %New interpolation
        eic = find(inpolygon(xc(:,1,i+1), xc(:,2, i+1), Xtri(e,:), Ytri(e,:))); %find which cells are in the element
        xce(eic) = e;
    end    
    c_xce = find(xce ~= 0);
    phi_c = zeros(length(c_xce),3);
    
    for e = 1:length(c_xce) %Interpolation
        phi_c(e,:) = [1, xc(c_xce(e) ,:, i+1)]*phi(:,:,xce(c_xce(e)));
    end
    
    X2C = sparse([c_xce, c_xce, c_xce], tr(xce(c_xce), :), phi_c, nc, np);
    
    PDGF_c = X2C*PDGF(:, i+1);
    fi_c = X2C*fi(:, i+1);
    
    %% Macrophage binding
    nmPDGF(:,i+1) = ( (1+ dt*n_bc)*eye(nm) + dt*n_gc*diag(PDGF_c(1:nm)))\(nmPDGF(:,i) + dt*n_gc*PDGF_c(1:nm));
    nmPDGF(death_m, i+1) = 0;
        
    %% TGF-concentration and TGF-receptors on fibroblasts  
    %Source function for TGF production
    sTGF = X2C(1:nm,:)'*(rtgf.*((1 - n_b)*fi_c(1:nm) + n_b));%*((n_a - 1)*nmPDGF(:, i+1) + 1)
    
    TGF(:,i+1) = (M + dt*(DTGF.max*Sc - (DTGF.max - DTGF.min)*F))\(M_old*TGF(:,i) + dt*sTGF );
    
    TGF_c = X2C(nm+1:end,:)*TGF(:,i+1);
    
    nftgf(:,i+1) = ((1+ dt*n_bc)*eye(nf) + dt*n_gc*diag(TGF_c))\(nftgf(:,i) + dt*n_gc*TGF_c);
    nftgf(death_f(nm+1:end),i+1) = 0;
    
    %% Collagen
    vfnc = zeros(nf, 2);
    vfnc(max(abs(vfc), [], 2) > 0, :) = diag(1./sqrt(sum(vfc(max(abs(vfc), [], 2) > 0, :).^2, 2))) * vfc(max(abs(vfc), [], 2) > 0, :); %Normalize velocity vectors for collagen synthesis
    
    colla14 = [collat(:,1,i); collat(:,4,i)]; %Save values of diagonals
    
    X2C_f = blkdiag(X2C(nm+1:end, :), X2C(nm+1:end, :)); %Interpolation matrix for coupled vector
    M_14 = blkdiag(M, M); %Mass matrix for coupled vector
    M_14_old = blkdiag(M_old, M_old);
    
    ru = [(rfc1 + (rfc2 - rfc1)*nftgf(:, i+1)).*vfnc(:,1).^2; (rfc1 + (rfc2 - rfc1)*nftgf(:, i+1)).*vfnc(:,2).^2]; %Receptor binding
    scolla14 = X2C_f'*spdiags(ru, 0, 2*nf, 2*nf); %Production of collagen to mesh points
    
    colla14 = (M_14 + dt*scolla14*X2C_f*[speye(np), speye(np); speye(np), speye(np)])\(M_14_old*colla14 + dt*scolla14*(1 - X2C_f*[fi(:, i+1); fi(:, i+1)])); %Calculate next iterate
    
    collat(:,1,i+1) = colla14(1:np); %Save in corresponding vectors
    collat(:,4,i+1) = colla14(np+1:end);
    collas(:,i+1) = collat(:,:,i+1)*[1;0;0;1];
    
    scolla2 = spdiags((rfc1 + (rfc2 - rfc1)*nftgf(:, i+1)).*(vfnc(:,1).*vfnc(:,2)), 0, nf, nf); %Calculate next iterate for remaining tensor elements
    collat(:,2,i+1) = M\(M_old*collat(:,2,i) + dt*X2C(nm+1:end, :)'*scolla2*(1 - X2C(nm+1:end, :)*fi(:, i+1) - X2C(nm+1:end, :)* collas(:,i+1)) );
    collat(:,3,i+1) = collat(:,2,i+1);
    
    %% Save mass matrix and contraction
    M_old = M;
    u_old = u;
    
    %% Report progress
    msg = sprintf('Processed %d/%d', i, T-1);
    fprintf([reverseStr, msg]);
    reverseStr = repmat(sprintf('\b'), 1, length(msg));
    
end
fprintf(' done\n');
toc

%% Plot solution
figure(1)

for i = 1:1:T
    subplot(2, 3, 1);
    trisurf(tr, xp(:,1,i), xp(:,2,i), tPA(:,i))
    title('tPA');
%         shading interp
    axis([xb, yb, 0, 1])
    caxis([0 1])
    view(180, 0)
    subplot(2, 3, 2);
    trisurf(tr, xp(:,1,i), xp(:,2,i), fi(:,i))
    title(['Fibrin ', num2str(ti(i))])
%         shading interp
    axis([xb, yb, 0, 1])
    view(3)
    caxis([0 1])
    subplot(2, 3, 4);
    trisurf(tr, xp(:,1,i), xp(:,2,i), TGF(:,i))
    title('TGF-\beta')
%         shading interp
    axis([xb, yb, 0, max(max(TGF))])
    view(3)
    caxis([0 1])
    subplot(2, 3, 5);
    trisurf(tr, xp(:,1,i), xp(:,2,i), PDGF(:,i))
    title('PDGF')
%         shading interp
    axis([xb, yb, 0, max(PDGF(:,1))])
    view(3)
    caxis([0 1])
    subplot(2,3,3)
    trisurf(tr, xp(:,1,i), xp(:,2,i), collas(:,i))
    title('Collagen')
%         shading interp
    axis([-3.1, 3.1, -2.1, 2.1, 0, 3])
    view(2)
    caxis([0 1])
    subplot(2,3,6)
    trisurf(tr, xp(:,1,i), xp(:,2,i), collas(:,i) + fi(:,i))
    title('Total fibers')
%         shading interp
    axis([xb, yb, 0, 1.1])
    view(3)
    caxis([0 1])
    %         pause
    drawnow;
end

%% Plot cell movement
figure(2)
Fv = zeros(np,4,T);

for i = 1:10:T
    clf    
    for j = 1:np
        [Dv, Fe] = eig(reshape(collat(j,:,i),2,2));
        Fei = Dv*Fe;
        Fv(j,:,i) = Fei(:)';
    end
    Fv(:,:,i) = diag(collas(:,i))*Fv(:,:,i);
    quiver(xp(:,1,i), xp(:,2,i), Fv(:,1,i), Fv(:,2,i), .2, 'k-', 'ShowArrowHead', 'off');
    hold on
    quiver(xp(:,1,i), xp(:,2,i), -Fv(:,1,i), -Fv(:,2,i), .2, 'k-', 'ShowArrowHead', 'off');    
    quiver(xp(:,1,i), xp(:,2,i), Fv(:,3,i), Fv(:,4,i), .2, 'k-', 'ShowArrowHead', 'off');  
    quiver(xp(:,1,i), xp(:,2,i), -Fv(:,3,i), -Fv(:,4,i), .2, 'k-', 'ShowArrowHead', 'off');  
    plot(xw(:,1,i), xw(:,2,i), 'r');
    plot(xp(gn,1,i), xp(gn,2,i), 'b');
    hold off
    
    axis([-3, 3, -2, 2.5], 'manual', 'equal')
    
    viscircles(xc(1:nm,:,i), rc*ones(nm,1), 'LineWidth', 1);
    viscircles(xc(nm+1:nc,:,i), rc*ones(nf,1), 'LineWidth', 1, 'EdgeColor', 'b');
    title(['Cell movement ', num2str(ti(i)), 'h'])
    drawnow;
end

%% Plot wound ratio
woundplotter = reshape(xw(end,:,:) - xw(1,:,:), 2, T)';
wounddist = sqrt(woundplotter(:,1).^2 + woundplotter(:,2).^2)/norm(xw(end,:,1) - xw(1,:,1));

xwn = [find(x(gn, 1) > xw(1,1), 1), find(x(gn, 1) < xw(end,1), 1, 'last')];
area = zeros(T,1);
for i = 1:T
    xwa = [xw(:,:,i); xp(gn(xwn(2):-1:xwn(1)), :, i)];
    area(i) = polyarea(xwa(:,1), xwa(:,2));
end
figure(5)
area = area/(area(1));
plot(ti, area, ti, wounddist);
legend('Area', 'Surface distance');

%% Label points
% figure(4)
% scatter(xp(:,1,end), xp(:,2,end), '.')
% labels = num2str((1:size(xp,1))','%d');
% axis([-3.1, 3.1, -2, 2.5], 'manual', 'equal')
% text(xp(:,1,end), xp(:,2,end), labels, 'horizontal','left', 'vertical','bottom')

%% Let computer sleep
% system('rundll32.exe powrprof.dll, SetSuspendState 0,1,0')