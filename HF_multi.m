% Script for solving the HF equations by energy minimization for the
% four-flavor model with a contact interaction

%% Parameters
clear param;

param.nf = 4;             % number of flavors
param.W  = 1;             % bandwidth
param.U  = 0.85;          % interaction strength
param.N =   2000;
param.E0 = 0.7;           % DOS parameters 
param.Evh = 0.7;
param.Ecut = 1/param.N;
param.nrand = 5;          % number of random initial guesses
param.B = 0;              % Zeeman field (deactivated)
param.spin = [1,-1,1,-1]';
param.is_lin_dos = 'y';     % if param.is_lin_dos=='y', use analytic expressions for linear DOS


%% Initialize DOS(eps), Ek(eps), n(eps)
eps = -param.W:param.W/param.N:param.W;
if (param.is_lin_dos == 'y')   % if not using linear DOS 
    % define DOS also for case is_lin_dos=='y' for compressibility computation
    nu = abs(eps)/param.W^2;
else        
    % Different DOS models
    %nu = abs(eps)/param.W^2;
    %nu = 0.1+abs(eps)/param.W^2.*sqrt((abs(eps)-param.Evh).^2+param.Ecut^2).^(-0.25);
    %nu = abs(eps)/param.W^2.*sqrt((abs(eps)-param.Evh).^2+param.Ecut^2).^(-0.25);
    %nu = abs(eps)/param.W^2;
    %nu = 0.1+abs(eps)/param.W^2.*log(sqrt((abs(eps)-param.Evh).^2+param.Ecut^2).^(-1));
    %nu = 0.1+abs(eps)/param.W^2;
    nu = abs(eps)/param.W^2.*(abs(eps)<=param.E0) ...
        + (param.W-abs(eps))/(param.W-param.E0)*param.E0/param.W^2.*(abs(eps)> param.E0);    
end
% Compute n as integral of nu(eps)
param.n = cumtrapz(eps, nu);
param.n = param.n - param.n(param.N+1);
param.n = param.n/param.n(end);

% Normalize nu such that n(eps=W)=1
ntot = trapz(eps, nu);
nu = nu*2/ntot;

% Compute Ek as integral of eps*nu
param.Ek = cumtrapz(eps, eps.*nu);
param.Ek = param.Ek - param.Ek(param.N+1);

[I,J] = meshgrid(1:param.nf, 1:param.nf);
param.utmat = I>J;


%% Loop over mu (chemical potential)
mu = 0:0.02:5;
%U   = 0.6:0.005:0.7;

n_min = zeros(param.nf, length(mu));
E_min = zeros(1, length(mu));
lb = [-1,-1,-1,-1]';
ub = [ 1, 1, 1, 1]';

n0 = [1,-1,0,0]';

options = optimoptions('fmincon');
options.Display = 'off';
options.OptimalityTolerance=1e-7;

Etr =   zeros(1,param.nrand+1);
ntr =  zeros(param.nf,param.nrand+1);
tic
for j = 1:length(mu)
    param.mu = mu(j);    % chemical potential for this iteration
    if (param.is_lin_dos == 'y')
        f = @(x)En1n2_lin(x,param); 
    else
        f = @(x)En1n2(x,param);  % Define f = En1n2 as anonymous function
    end
    [ntr(:,1),Etr(1)] = fmincon(f,n0,[],[],[],[],lb,ub,[],options);
    % Random initial conditions
    n0tr = 2*(rand(param.nf,param.nrand)-0.5);
    for jr=1:(param.nrand+1)        
        if (jr==1)
            n_init = n0;
        else
            n_init = n0tr(:,jr-1);
        end
        [ntr(:,jr),Etr(jr)] = fmincon(f,n_init,[],[],[],[],lb,ub,[],options);
    end
    [~,I] = min(Etr); 
    n_min(:,j) = ntr(:,I);
    E_min(j) = Etr(I);
    n0 = n_min(:,j);
    if (mod(j,50)==0)
        fprintf('mu = %.4f\n',mu(j))              
    end
end

n_min = sort(n_min,1);
n_tot = sum(n_min,1);

%% Compute compressibility
nu_a = interp1(param.n,nu,n_min(:)).*(abs(n_min(:))<0.999);
nu_a = reshape(nu_a, size(n_min,1), size(n_min,2));
nu_bar = sum(nu_a./(1-param.U*nu_a),1);
dndm = nu_bar./(1 + param.U*nu_bar);



%% Plot the results
figure
subplot(2,2,1);
plot(eps,nu,'linewidth',1);
prettyfig; xlabel('$$\varepsilon/W$$'); ylabel('$$\nu$$');
subplot(2,2,2);
plot(mu, n_min,'linewidth',1)
prettyfig; xlabel('$$\mu$$'); ylabel('$$n_\alpha$$');
subplot(2,2,3);
plot(n_tot,mu,'linewidth',1); 
prettyfig; xlabel('$$n$$'); ylabel('$$\mu$$');
subplot(2,2,4);
%plot(n_tot(2:end), diff(mu)./diff(n_tot),'linewidth',1)
plot(n_tot, 1./dndm,'linewidth',1)
prettyfig; xlabel('$$n$$'); ylabel('$$d \mu/dn$$');
axis([0 4 0 3])

toc







