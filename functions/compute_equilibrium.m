function [equilibrium_outcomes, simulations_outcomes, agg_stat] = compute_equilibrium(param_struct)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Project: Risk-sharing in a dual market
% Créchet (2020)
% matlab function
% file name: "compute_equilibrium.m"
% created: 28-10-2020
% last updated: sept 2023
% Computes numerical equilibrium solution

% argument: a structure with (i) model economic parameters; (ii) parameters
% specifying options for numerical algorithm; (iii)  options for
% simulations.

% out: equilibrium_outcomes: equilibrium variables (e.g., value and policy functions)
% simulation_outcomes: panel with simulated individual employment histories
% agg_stat: equilibrium statistics

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Parameters

% timer start
timer = tic; % timer start

% stack out parameters from struture param
p = param_struct;
ind = p.ind;
vfi_tol = p.vfi_tol;
vfi_tol_U = p.vfi_tol_U;
equilibrium = p.equilibrium;
compute_welfare = p.compute_welfare;
estimate_density = p.estimate_density;

pval = p.pval;
Ks = p.sim.Ks;
Iz = p.Iz;
Ix = p.Ix;
rng_algo = p.sim.rng_algo;
max_iter = p.max_iter;

% economic parameters
% preferences
bbeta = pval(ind.beta);
ssigma = pval(ind.sigma);

% matching
eeta = pval(ind.eta);
A = pval(ind.A);
kkappa = pval(ind.kappa);

% prod. stoch processes
sigma_x = pval(ind.sigma_x);
llambda = pval(ind.lambda);
zub = pval(ind.zub);
zlb = pval(ind.zlb);
ddelta = pval(ind.delta);

% employment protection/u benef.
F = pval(ind.F);
pphi = pval(ind.phi);
pphi0 = pval(ind.phi0);
b = pval(ind.b);

% reservation wage (a param. if partial equil.)
wr = pval(ind.wr);

% LM tightness (normalize to one if partial equil.)
theta = pval(ind.theta);

% options for fzero
opts = optimset('Display', 'notify');

% (note: partial equil: b is set to be consistent with reservation wage and
% tightness)

%% 0a. Preferences and production technology: functional forms

% production function
y = @(z,x) (z .* x);

% utility function, inverse and derivative
if ( ssigma ~= 1 )
    u =       @(c) (c.^(1-ssigma)-1)./(1-ssigma);
    u_inv =   @(u) (1-(ssigma-1).*u).^(1/(1-ssigma));
else
    u =       @(c) log(c);
    u_inv =   @(u) exp(u);
end

if ssigma == 0
    u = @(c) c;
    u_inv = @(u) u;
end

u_prim =  @(c) c.^(-ssigma);

%% 0b. Grids for exogenous variables

% grid and pmf for stochastic output component
zgrid = linspace(zlb, zub, Iz)';
gz = (1/Iz)*ones(Iz,1);

% grid and distribution for the match quality
mu_x = - sigma_x^2/2; % unconditional mean equals one
[xtmp, gx] = n_discretized(mu_x, sigma_x, Ix, 0.99);
xgrid = exp(xtmp);


%% 1a. Compute maximum wages in an optimal PC

% maximum continuation wages
tic
wz = zeros(Iz,Ix);

for ix = 1:Ix

    % match quality
    x = xgrid(ix);

    for iz = 1:Iz

        % stochastic shock
        z = zgrid(iz);

        % compute probability of a equal or lower prod cond. on a shock
        Pr = gz'*(zgrid<=z);

        % expected match output cond. on z'>z
        if iz < Iz
            Ez = gz'*(y(zgrid,x).*(zgrid>z)) / (1-Pr);
        else
            Ez = 0;
        end

        % present value of expected lifetime output net of expected firing costs
        Y = (1-Pr)/(1-bbeta*(1-ddelta)*(1-llambda*Pr)) * (Ez - bbeta*(1-ddelta)*llambda*Pr*F);

        % compute max wage
        K = 1 + bbeta*(1-ddelta)*llambda*(1-Pr)/(1-bbeta*(1-ddelta)*(1-llambda*Pr));
        wz(iz,ix) = K^(-1)*( y(z,x) + bbeta*(1-ddelta)*llambda*(Y - Pr*F) ...
            + (1-bbeta*(1-ddelta)*(1-llambda))*F );

    end
end

% compute PC profit function cond. on match quality
Jp = cell(Ix,1);
for ix = 1:Ix

    % match quality
    x = xgrid(ix);

    % probability of a wage cut/PC binding
    Pr = @(w) gz'*(w > wz(:,ix));

    % expected discounted output
    EY = @(w,z) (1-bbeta*(1-ddelta)*(1-llambda))^(-1)*( y(z,x) + ....
        + bbeta*(1-ddelta)*llambda / (1-bbeta*(1-ddelta)*(1-llambda*Pr(w)))*( (gz.*(w <= wz(:,ix)))'*y(zgrid,x)) );

    % profit function
    Jp{ix} = @(w,z) EY(w,z) - (1 - bbeta*(1-ddelta)*(1-llambda*Pr(w)))^(-1)*( w + bbeta*(1-ddelta)*llambda*Pr(w)*F );

end

% % check if wz and Jp are consistent (should have Jp(wz(z,x),z,x) = -F)
% g = @(iz,ix) Jp{ix}(wz(iz,ix), zgrid(iz)) + F;
% G = zeros(Iz, Ix);
% for ix = 1:Ix
%     for iz = 1:Iz
%         G(iz,ix) = g(iz,ix);
%     end
% end
% disp(max(abs(G), [], 'all'))


% solve for max wage in a starting PC using Jp function
w0z = wz;
if F>0
    for ix = 1:Ix
        for iz = 1:Iz
            fn = @(w) Jp{ix}(w,zgrid(iz));
            w0z(iz,ix) = fzero(fn, wz(iz,ix), opts);
        end
    end
end

% % could check that w0z = wz when F is set to F=eps
% disp( max( abs(wz-w0z), [], 'all' ) )


%% 1b. Compute maximum wage of an optimal temp contract
if F>0

    wtz = zeros(Iz, Ix);

    %parfor ix = 1:Ix
    parfor ix = 1:Ix

        x = xgrid(ix);

        for iz = 1:Iz

            % stochastic shock
            zz = zgrid;
            z = zz(iz);

            % probability of equal or lower prod cond. on a shock
            Pr = gz'*(zgrid<=z);

            % expected match output conditional on z'>z
            if iz < Iz
                Ez = (gz.*(zgrid>z))'*y(zgrid,x) / (1-Pr);
            else
                Ez = 0;
            end

            % expected PC conversion profit
            EJp_u = @(w) gz'*max(Jp{ix}(w*ones(Iz,1), zgrid), 0);

            % expected PC profits, conditional on z'>z
            if iz<Iz
                EJp_c = @(w) (gz.*(zgrid>z))'*max(Jp{ix}(w*ones(Iz,1), zgrid), 0) / (1-Pr);
            else
                EJp_c = @(w) 0;
            end

            % function of (w,z) to find root for
            K = bbeta*(1-ddelta)*(1-pphi)*llambda*(1-Pr) / (1-bbeta*(1-ddelta)*(1-pphi)*(1-llambda));

            g = @(w) y(z,x) + K/(1-K)*Ez - 1/(1-K)*w ...
                + bbeta*(1-ddelta)*pphi*( (1-llambda)*max(Jp{ix}(w,z),0) + llambda*EJp_u(w) ) ...
                + K/(1-K)*bbeta*(1-ddelta)*pphi*( (1-llambda)*EJp_c(w) + llambda*EJp_u(w)) ;

            % solve for root and deduce max wage for (z,x)
            lb = w0z(iz, ix); ub = wz(iz, ix);
            if abs(g(lb)) < 1e-10
                wtz(iz,ix) = lb;
            else
                wtz(iz,ix) = fzero(g, [lb ub], opts);
            end

        end
    end
end
time.max_wages = toc;

% % could check that wtz = wz when F is set to F=eps
% disp( max( abs(wtz-wz), [], 'all' ) )



%% 2.a. Initialize equilibrium loop

equil_timer = tic;

% grid for max current utility
uz = u(wz);
u0z = u(w0z);
if F>0
    utz = u(wtz);
end

% reservation utility
ur = u(wr);

% Value of unemployment consistent with reservation wage
U = ur/(1-bbeta);

% (note: if gen. equilibrium: U is the initial guess on equilibrium U
% lifetime value)

% initialize VFI algo (if gen. equilibrium)
U0 = U;
dist_U = inf;
iter = 0;


%% loop for VFI iterations on U
while (dist_U > vfi_tol_U && iter <= max_iter )


    %% 2b. Compute maximum LTU in optimal permanent contract
    tic

    % Compute max continuation LTU values with VFI
    Vz = ones(Iz, Ix)*U;

    for ix = 1:Ix
        V0 = 0;
        dist = inf;

        while dist>vfi_tol

            EV = (1-llambda)*Vz(:,ix) + llambda*min(Vz(:,ix), Vz(:,ix)')*gz;
            Vz(:,ix) = uz(:,ix) + bbeta*(1-ddelta)*EV + bbeta*ddelta*U;
            Vz(:,ix) = max(Vz(:,ix), U);
            dist = max(abs(Vz(:,ix)-V0), [], "all");
            V0 = Vz(:,ix);

        end
    end

    % Deduce worker's LTU as function of wage and match state
    Vp = cell(Ix,1);
    for ix = 1:Ix
        Pr = @(u) gz'*(u > uz(:,ix));
        Vp{ix} = @(u) max( (1-bbeta*(1-ddelta)*(1-llambda*Pr(u)))^(-1)*( u ...
            + bbeta*(1-ddelta)*llambda*( (gz.*(u > uz(:,ix)))'*Vz(:,ix) ) + bbeta*ddelta*U ), U);
    end

    %     % check consistency
    %     g = @(iz,ix) Vp{ix}(uz(iz,ix)) - Vz(iz,ix);
    %     G = zeros(Iz,Ix);
    %     for ix = 1:Ix
    %         for iz = 1:Iz
    %         G(iz,ix) = g(iz,ix);
    %         end
    %     end
    %     disp(max( abs(G), [], 'all'))


    % deduce max LTU for a starting PC
    V0z = Vz;
    if F>0
        for ix = 1:Ix
            for iz = 1:Iz
                V0z(iz,ix) = max(Vp{ix}(u0z(iz,ix)), U);
            end
        end
    end


    %% 3b. Max LTU in optimal temp. contract
    if F>0

        % max LTU
        Vtz = ones(Iz,Ix)*U;

        for ix = 1:Ix

            V0 = 0;
            dist = inf;

            while dist>vfi_tol

                % expectation cond. on continuation in temp. contract
                EVt = (1-llambda)*Vtz(:, ix) + llambda*min( Vtz(:,ix), Vtz(:,ix)')*gz;

                % conditional on conversion in perm.
                vp = min( Vp{ix}(utz(:,ix)),  V0z(:,ix) );
                tmp =  min( vp,  V0z(:,ix)' );
                EVp = (1-llambda)*vp + llambda*(tmp*gz);

                % unconditional
                EV = (1-pphi)*EVt + pphi*EVp;

                % VF
                Vtz(:,ix) = utz(:,ix) + bbeta*(1-ddelta)*EV + bbeta*ddelta*U;
                Vtz(:,ix) = max(Vtz(:,ix), U);
                dist = max(abs(Vtz(:,ix)-V0), [], "all");
                V0 = Vtz(:,ix);

            end
        end

    end

    %     % check
    %     disp( max( abs(Vtz - Vz), [], 'all' ) )
    %     time.max_ltu = toc;


    %% 4. Matching surplus and hiring policy

    tic
    % wage grid
    if F>0
        tmp = sort([wz(Iz,:), w0z(Iz,:), wtz(Iz,:)])';
    else
        tmp = wz(Iz,:)';
    end
    wgrid = unique( [wr; tmp(tmp>wr)] );
    ugrid = u(wgrid);
    Iv = length(wgrid);


    % Grids for PC functions (for welfare analysis or for computing TC surplus)
    if F>0 || compute_welfare

        % Grid for PC profits on a grid

        Jp_grid = zeros(Iv, Iz, Ix);
        for ix = 1:Ix
            for iv = 1:Iv
                w = wgrid(iv);
                Jp_grid(iv,:,ix) = Jp{ix}(w,zgrid);
            end
        end

        % worker's surplus, VFI
        Vp_grid = zeros(Iv, Iz, Ix);
        for ix = 1:Ix
            for iv = 1:Iv
                uu = ugrid(iv);
                Vp_grid(iv,:,ix) = min( Vp{ix}(uu), V0z(Iz,ix) );
            end
        end
    end


    %% 4a. PC matching surplus

    % employer
    Se_p = zeros(Iv,Ix);
    for ix = 1:Ix
        for iv = 1:Iv
            Se_p(iv, ix) = max( Jp{ix}(wgrid(iv), zgrid(Iz)), 0);
        end
    end

    % worker
    Sw_p = zeros(Iv,Ix);
    for ix = 1:Ix
        for iv = 1:Iv
            Sw_p(iv, ix) = max( min( Vp{ix}(ugrid(iv)), V0z(Iz, ix) ) - U, 0);
        end
    end

    % total surplus of a perm job given the wage and match quality (in numéraire units)
    Sp = Se_p + (u_prim(wgrid).^(-1)).*Sw_p;



    %% 4b. TC matching surplus

    if F>0

        % Employer (VFI)
        % Init array for TC profits grid
        Jt_grid = Jp_grid;

        % loop over match quality
        parfor ix = 1:Ix

            % initialize
            dist = inf;
            Jt0 = 0;

            % current profit grid
            pi_grid = y(zgrid, xgrid(ix))' - wgrid;

            % VFI
            while dist > vfi_tol

                % next period expectation
                tmp = (1-pphi)*max(Jt_grid(:,:,ix), 0) + pphi*max(Jp_grid(:,:,ix), 0);
                EJt = (1-llambda)*tmp + llambda*tmp*gz;
                Jt_grid(:,:,ix) = pi_grid + bbeta*(1-ddelta)*EJt;

                % update
                dist = max( abs(Jt_grid(:,:,ix)-Jt0), [], "all" );
                Jt0 = Jt_grid(:,:,ix);

            end
        end

        % Deduce matching surplus
        Se_t = reshape( max(Jt_grid(:,Iz,:), 0), [Iv Ix] );

        % % check (when F=eps)
        % disp(max( abs(Se_t - Se_p), [], "all" ))


        % VFI
        Vt_grid = Vp_grid;

        parfor ix = 1:Ix

            dist = inf;
            V0 = 0;

            while dist>vfi_tol

                % expectation
                tmp = (1-pphi)*min(Vt_grid(:,:,ix), Vtz(:,ix)') + pphi*min(Vp_grid(:,:,ix), V0z(:,ix)');
                EV = (1-llambda)*tmp + llambda*tmp*gz;

                % updated VF
                Vt_grid(:,:,ix) = ugrid + bbeta*(1-ddelta)*EV + bbeta*ddelta*U;

                % convergence
                dist = max(abs(Vt_grid(:,:,ix)-V0), [], "all");
                V0 = Vt_grid(:,:,ix);

            end
        end

        % matching surplus
        tmp = reshape( Vt_grid(:,Iz,:), [Iv Ix] );
        Sw_t = max( min(tmp, Vtz(Iz,:)) - U, 0);

        % total (in numéraire units)
        St = Se_t + (u_prim(wgrid).^(-1)).*Sw_t;

        % % check
        % disp(max( abs(St - Sp), [], "all" ))

    end


    %% 4c. Optimal hiring policy

    % hiring cutoff, PC
    ix = max(sum( w0z(Iz, :) < wr ), 1);
    if ix < Ix
        tmp = griddedInterpolant( xgrid, w0z(Iz, :), 'linear');
        fn = @(x) tmp(x) - wr;
        xp = fzero( fn, xgrid(ix), opts);
    else
        xp = inf;
    end


    % TC
    if F>0
        ix = max(sum( wtz(Iz, :) < wr ), 1);
        if ix < Ix
            tmp = griddedInterpolant( xgrid,  wtz(Iz, :), 'linear');
            fn = @(x) tmp(x) - wr;
            xt = fzero( fn, xgrid(ix), opts);
        else
            xt = inf;
        end
    end

    % % check
    % disp(xt-xp)

    % Expected surplus
    if F>0

        % grid matrix for hiring policy decision
        D = Sp - St;
        D( abs(D) < 1e-8 ) = 0; % assign zero when diff is approx zero to some tol.
        Hp = (1-pphi0)*(D>=0) + pphi0*(Sp>0);
        Ht = (1-pphi0)*(D<0 & St>0);

        % expected surplus grid
        ESe_grid = (Hp.*Se_p + Ht.*Se_t)*gx;
        ESw_grid = (Hp.*Sw_p + Ht.*Sw_t)*gx;

    else

        ESe_grid = Se_p*gx;
        ESw_grid = Sw_p*gx;

    end

    % Interpolate expected surplus functions
    ESe = griddedInterpolant(wgrid, ESe_grid, 'linear');
    ESw = griddedInterpolant(wgrid, ESw_grid, 'linear');

    % FOC for bargained wage
    g = @(w) (1-eeta)*(u_prim(w).^(-1)).*ESw(w) - eeta.*ESe(w);

    % solve for bargained wage
    iw = sum( g(wgrid) < 0 );
    if (iw > 0 && iw < length(wgrid))
        lb = wgrid(iw);
        ub = wgrid(iw+1);
        w_new = fzero(g, [lb ub]);
    else
        w0 = wr;
        w_new = fzero(g, w0, opts);
    end


    %     % check: find wage that max weighted surplus; it should be equivalent
    %     to the one found using the FOC.

    %     tmp = ESe_grid.^(1-eeta).*ESw_grid.^eeta;
    %     Omega = griddedInterpolant(wgrid, tmp, 'spline');
    %     [~, iw] = max(tmp);
    %     fn = @(w) -Omega(w);
    %     opts_fmincon = optimoptions('fmincon', 'Display', 'notify');
    %     w_new2 = fmincon( fn, wgrid(iw), wgrid(1), wgrid(end), [], [], [], [], [], opts_fmincon);


    % deduce agents' initial matching values
    Jnew = ESe(w_new);
    Vnew = ESw(w_new) + U;


    % Optimal PC/TC choice
    if F>0

        % interpolate marginal benefit function
        [ww, xx] = ndgrid(wgrid, xgrid);
        tmp = griddedInterpolant(ww, xx, D, 'linear');

        % 'inaction' cutoff
        ix = sum( 1 - ( D(1,:) >= 0 & xgrid' >= xp ) ) + 1;
        if ix <= Ix
            xtilde = xgrid(ix);
        else
            xtilde = nan;
        end

        % marginal benefit as a function of x, and given bargained wage
        fn = @(x) tmp(w_new, x);

        % case 1: inaction cutoff xtilde exists
        if ~isnan(xtilde)

            % case 1a: xtilde is the solution (e.g., when worker BP = 0)
            if fn(xtilde)<1e-8

                xhat = xtilde;

                % case 1b: xhat < xtilde
            else
                xhat = fzero(fn, [xp xtilde], opts);
            end

        end

        % case 2: inaction cutoff xtilde does not exist
        if isnan(xtilde)

            % 2a: choice cutoff xhat exists
            if fn(xp)*fn(xgrid(Ix)) < 0
                xhat = fzero(fn, [xp xgrid(Ix)], opts);

                % 2b: choice cutoff does not exist
            else
                xhat = inf;
            end

        end
    end


    % compute the expected mean hiring wage

    % wage cond. on match quality and perm
    wp = min( w_new, w0z(Iz,:) )';

    if F>0

        % hiring policy vectors
        hp = (1-pphi)*(xgrid>=xhat) + pphi*(xgrid>=xp);
        ht = (1-pphi)*(xgrid>=xt & xgrid<xhat);

        % conditional wage in temp
        wt = min( w_new, wtz(Iz,:) )';

        % probability of hiring
        P = (hp + ht)'*gx;

        % distribution conditional on hiring
        gx_new = (gx.*(hp + ht))/P;

        % average log wage conditional on hiring
        wnew_mn = (hp.*log(wp) + ht.*log(wt))' * gx_new;

    else

        hp = (xgrid>=xp);
        P = hp'*gx;
        gx_new = (gx.*hp)/P;
        wnew_mn = log(wp)'*gx_new;

    end

    clearvars wp wt
    time.hiring_pol = toc;


    %% 5. Update unemployment value (gen. equil) or backout b and kappa (partial equil.)

    % partial equilibrium: compute kappa and b consistent with preset values for wr and theta
    if strcmp(equilibrium,'partial')

        % back out vacancy posting cost consistent with theta = 1
        kkappa = bbeta*A*theta^(-eeta)*Jnew;

        % back out non-work income consistent with reservation wage and
        % theta = 1
        b = u_inv( (1-bbeta)*U - bbeta*A*theta^(1-eeta)*(Vnew-U) );

        % no additional iteration: set crit_u = 0
        dist_U = 0;

    end

    % general equilibrium: update unemployment LTU, U
    if strcmp(equilibrium,'general')

        iter = iter + 1;

        % compute labor market tightness consistent with profits
        theta = (bbeta*A*Jnew/kkappa)^(1/eeta);

        % update value of unemployment given the tightness
        U = u(b) + bbeta*(A*theta^(1-eeta)*(Vnew-U0) + U0);

        % update reservation wage
        wr = u_inv( (1-bbeta)*U );

        % check convergence
        dist_U = norm(U-U0);
        U0 = U;
        disp(dist_U);

        % display a message when convergence do occur given tol.
        if dist_U < vfi_tol_U
            disp('U converges');
        end

    end


end % end of loop on U iterations

% notify if no convergence given tol.
if iter >= max_iter
    disp('no convergence of U given tol, dist:')
    disp(dist_U)
end

% save equilibrium outcomes
eql = ...
    struct( 'U', U, 'wr', wr, 'theta', theta, 'Jnew', Jnew, 'Vnew', Vnew, 'Vz', Vz, ...
    'w_new', w_new, 'wz', wz, 'xp', xp, 'b', b, 'kkappa', kkappa, 'dist_U', dist_U, ...
    'xgrid', xgrid, 'zgrid', zgrid, 'Sp', Sp, 'wgrid', wgrid, 'ESe', ESe, 'ESw', ESw);

if F>0
    eql.Vtz = Vtz;
    eql.V0z = V0z;
    eql.D = D;
    eql.xhat = xhat;
    eql.xtilde = xtilde;
    eql.xt = xt;
    eql.wtz = wtz;
    eql.w0z = w0z;
    eql.St = St;
    eql.D = D;
end

time.equilibrium = toc(equil_timer);

%% 6. Simulations

tic

% wage law of motions
% interpolate maximum wage policies
[zz, xx] = ndgrid(zgrid,xgrid);
w0z_grid = w0z;
wz_grid = wz;
w0z = griddedInterpolant(zz, xx, w0z_grid, 'linear');
wz  = griddedInterpolant(zz, xx, wz_grid, 'linear');

% permanent contract, continuation
wp = @(w,z,x) min(w, wz(z,x));

% permanent, hiring/conversion
wp0 = @(w,z,x) min(w, w0z(z,x));

% temporary contract
if F>0
    wtz_grid = wtz;
    wtz  = griddedInterpolant(zz, xx, wtz_grid, 'linear');
    wt = @(w,z,x) min(w,wtz(z,x));
else
    wt = [];
end

% reservation wage
wr = u_inv( (1-bbeta)*U );

% separation rule
dp  = @(z,x) (wz(z,x) < wr);
dp0 = @(z,x) (w0z(z,x) < wr);
if F>0
    dt  = @(z,x) (wtz(z,x) < wr);
else
    dt = [];
end

% simulate employment histories
% draw random variables
% match quality: log normal with par. mu_x, sigma_x
% initialize stream
s0 = RandStream(rng_algo, "Seed", 0);
xsim = exp( sigma_x*randn(s0, [Ks,1]) + mu_x );

% ind_tmp = discretize( xsim, xgrid, 'IncludedEdge', 'right' );
% xsim_discr = xgrid( ind_tmp );

% 'hiring-regulation' shock
s0 = RandStream(rng_algo, "Seed", 1);
if F>0
    no_TC = (rand(s0, [Ks,1])<=pphi0);
else
    no_TC = ones(Ks,1);
end

% dummy for hiring decisions
if F==0
    hiring = (xsim>=xp);
    hiring_perm = hiring;
else
    hiring = (no_TC & (xsim>=xp)) | (no_TC==0 & (xsim>=xt));
    hiring_perm = (no_TC & (xsim>=xp)) | (no_TC==0 & (xsim>=xhat));
end

% keep only vector entries for match quality leading to an employemnt spell
xsim = xsim( hiring );
hiring_perm = hiring_perm( hiring );

% number of employment spells
Lx = length( xsim );

% draw max duration of emp. spells based on exog. separation prob.
s0 = RandStream(rng_algo, "Seed", 2);
dur = @(u) log(1 - u)./(log(1 - ddelta)) - 1;
tmp = rand(s0, [1 Lx]);
lz = int64( max( ceil(dur(tmp)), 1) );

% draw stochastic shocks
s0 = RandStream(rng_algo, "Seed", 3);
Lz = max(lz);
zsim = zeros(Lz,Lx);
zsim(1,:) = 1;

for ix = 1:Lx
    shock = ( rand(s0, [lz(ix) 1]) <= llambda );
    zz = zgrid( randi(s0, [1, Iz], [lz(ix) 1]) );

    for iz = 2:lz(ix)
        zsim(iz,ix) = (1-shock(iz))*zsim(iz-1,ix) + shock(iz)*zz(iz);
    end
end

% initialize matrices for wage, employment status, tenure
emp  = false(Lz,Lx);
perm = false(Lz,Lx);
wage = zeros(Lz,Lx);
ten  = int64(zeros(Lz,Lx));

% initial employment status
emp(1,:) = true;
perm(1,:) = hiring_perm';

%  hiring wage
xx = xsim; zz = ones(Lx,1);
if F==0
    wmax = w0z(zz,xx)';
else
    wmax = perm(1,:).*w0z(zz,xx)' + (1-perm(1,:)).*wtz(zz,xx)';
end
wage0 = min( w_new, wmax );
wage(1,:) = wage0;

% loop over employment spells
parfor ix = 1:Lx

    % match quality
    x = xsim(ix);

    % initialize predetermined variables (state + tenure)
    w0 = wage0(ix);
    perm0 = hiring_perm(ix);
    ten0 = int64(0);

    % random stream for employment spells
    s1 = RandStream(rng_algo, "Seed", ix+3);

    for iz = 2:Lz

        % stochastic component
        z = zsim(iz,ix);

        % case 1: permanent contract
        if (perm0==true && z>0)

            emp(iz,ix) = (1-dp(z,x));
            perm(iz,ix) = emp(iz,ix);
            wage(iz,ix) = emp(iz,ix)*wp(w0, z, x);

        elseif (perm0==false && z>0)

            % case 2: temporary contract
            conversion = (rand(s1)<=pphi);
            perm(iz,ix) = conversion*(1-dp0(z,x));
            emp(iz,ix) = perm(iz,ix) + (1-conversion)*(1-dt(z,x));
            wage(iz,ix) = emp(iz,ix)*( perm(iz,ix)*wp0(w0, z, x) + (1-perm(iz,ix))*wt(w0, z, x) );

        end

        % update tenure
        ten(iz,ix) = int64(emp(iz,ix))*(ten0 + 1);

        % update pred. state variables
        w0 = wage(iz,ix);
        perm0 = perm(iz,ix);
        ten0 = ten(iz,ix);

        % exit employment spell loop after match separation
        if ( emp(iz,ix) == 0 )
            break
        end

    end
end

% stack in vectors, remove zeros (to keep only employment spell)
LL = Lx*Lz;
emp = reshape( emp, [LL 1] );

% match quality
xsim = repmat( xsim', [Lz 1]);
xsim = reshape( xsim, [LL 1]);
xsim = xsim(emp>0);

% output shock
zsim = reshape( zsim, [LL 1]);
zsim = zsim(emp>0);

% wage, tenure, contract type
wage = reshape( wage, [LL 1]);
wage = wage(emp>0);
ten  = reshape( ten, [LL 1]);
ten = ten(emp>0);
perm = reshape( perm, [LL 1]);
perm = perm(emp>0);

% adjust employment vector (to get the size of simulated emp. spell
% vectors)
% emp = emp(emp>0);

% match productivity
ysim = zsim.*xsim;

% table
simulations_outcomes = table(perm, zsim, xsim, ysim, wage, ten);
time.simulations = toc;



%% 6b. estimate match-quality distributions and policy functions
tic

L = 100;
xbins = exp( linspace(log(min(xsim)), log(max(xsim)), L) )';

if strcmp( estimate_density, 'kernel')

    hp = ksdensity( log(xsim(perm)), log(xbins), "BoundaryCorrection", "reflection", "Support", [log(xp), Inf]);

    if F>0 && pphi0<1
      
        % unconditional
        h = ksdensity( log(xsim), log(xbins), "BoundaryCorrection", "reflection", "Support", [log(xt), Inf]);

        % temp
        prob_perm = sum(perm)/sum(emp);
        ht = (h - hp*prob_perm) / (1-prob_perm);
        
    else
        ht = 0*hp;
        h = hp;
    end

    % domain of separation-probability function
    xbins_sep = xbins;

else

    % permanent
    np = histcounts(xsim(perm), xbins)';
    hp = np/sum(np);

    % temp
    if F>0

        % unconditional
        n = histcounts(xsim, xbins)';
        h = n/(sum(n));

        % temp
        prob_perm = sum(perm)/sum(emp);
        ht = (h - hp*prob_perm) / (1-prob_perm);

    else

        ht = 0*hp;
        h = hp;
       
    end

    % domain of separation function
    xbins_sep = zeros(L-1,1);
    for i = 1:L-1
        xbins_sep(i) = (xbins(i) + xbins(i+1))/2;
    end

end

% separation proba function
% perm
% prob of z'<zp cond. on x and on a shock
Gz = @(x) gz'*dp(zgrid,x*ones(Iz,1));
sp_f = @(x) ddelta + (1-ddelta)*llambda*Gz(x);

% separation probability over bins
sp = zeros(length(xbins_sep),1);
for i = 1:length(xbins_sep)
    if hp(i) > 0
        sp(i) = sp_f(xbins(i));
    end
end

% temporary: separation cond. on z,x
if F>0

    % prob of z'<zc cond. on x and on a shock
    Gz0 = @(x) gz'*dp0(zgrid,x*ones(Iz,1));

    % prob of z'<zt cond. on x and on a shock
    Gzt = @(x) gz'*dt(zgrid,x*ones(Iz,1));

    % function for fraction of matches with z'<z_c (cond. on x)
    H = @(x) (x<=xp) * 1 ...
        + (x > xp) * (1-ddelta)*(1-pphi)*llambda*( Gz0(x) - Gzt(x) ) /  ( ddelta + (1-ddelta)*( (1-pphi)*llambda*(1-Gzt(x)) + pphi ));

    % separation function
    st_f = @(x) ddelta + (1-ddelta)*(1-pphi)*llambda*((1-pphi)*Gzt(x) + pphi*Gz0(x)) + (1-ddelta)*(1-llambda)*pphi*H(x);

    % conversion probability function
    c_f = @(x) (1-ddelta)*pphi*( (1-llambda)*(1-H(x)) + llambda*(1-Gz0(x)) );

    % simulated
    st = zeros(length(xbins_sep),1);
    c = zeros(length(xbins_sep),1);
    for i = 1:length(xbins_sep)
        if ht(i) > 0
            st(i) = st_f(xbins(i));
            c(i) = c_f(xbins(i));
        end
    end
else
    st = zeros(length(xbins_sep),1);
    c = [];
end

eql.h = h;
eql.sp = sp;
eql.hp = hp;
if F>0
    eql.c = c;
    eql.st = st;
    eql.ht = ht;
end
if strcmp(estimate_density, 'kernel')
    eql.x = xbins;
else
    eql.x = xbins_sep;
end

time.distr_pol_func = toc;


%% 6c. Simulated welfare

if compute_welfare

    % Interpolate permanent job functons over state space
    [ww, zz, xx] = ndgrid(wgrid, zgrid, xgrid);
    uu = u(ww);
    Jp = griddedInterpolant(ww, zz, xx, Jp_grid, 'linear');
    Vp = griddedInterpolant(uu, zz, xx, Vp_grid, 'linear');

    % temporary jobs
    if F>0
        Jt = griddedInterpolant(ww, zz, xx, Jt_grid, 'linear');
        Vt = griddedInterpolant(uu, zz, xx, Vt_grid, 'linear');
    end

    % workers' expected lifetime utility value and firms profits
    util = u(wage);
    if F==0
        simulations_outcomes.V = Vp(util, zsim, xsim);
        simulations_outcomes.J = Jp(wage, zsim, xsim);
    else
        simulations_outcomes.V = perm.*Vp(util, zsim, xsim) + (1-perm).*Vt(util, zsim, xsim);
        simulations_outcomes.J = perm.*Jp(wage, zsim, xsim) + (1-perm).*Jt(wage, zsim, xsim);
    end

    % in lifetime consumption equivalent terms
    simulations_outcomes.cV = u_inv( (1-bbeta)*simulations_outcomes.V );

    % lifetime annuity equivalent
    simulations_outcomes.cJ = (1-bbeta)*simulations_outcomes.J;

end


%% 7. Final statistics
tic

% Job finding rate
G = @(x) logncdf(x, mu_x, sigma_x);
if F>0
    UE = A*theta^(1-eeta)*( (1-pphi0)*(1-G(xt)) + pphi0*(1-G(xp)) );
    UP = A*theta^(1-eeta)*( (1-pphi0)*(1-G(xhat)) + pphi0*(1-G(xp)) );
else
    UE = A*theta^(1-eeta)*(1-G(xp));
    UP = UE;
end
UT = UE - UP;

% separation rate
if strcmp(estimate_density, 'kernel')
    PU = trapz( log(xbins), hp.*sp );
else
    PU = hp'*sp;
end

if F>0
    if strcmp(estimate_density, 'kernel')
        TU = trapz( log(xbins), ht.*st );
        TP = trapz( log(xbins), ht.*c );
    else
        TU = ht'*st;
        TP = ht'*c;
    end
else
    TU = nan;
    TP = nan;
end

% temp employment share and unconditional separation rate
if F>0 && pphi0<1
    T = 1-prob_perm;
    %T = UT / (TU + TP) * (U/E);
    EU = (1-T)*PU + T*TU;
else
    EU = PU;
    T = 0;
end

% unemployment rate
U = EU / ( EU + UE );
E = 1-U;
%U = PU*(TU+TP) / ( (TU+TP)*(PU+UT+UP) - UT*(TU-PU) );

% Temp employment flows
T_inflow = UT/UE;
T_outflow = T*TU/EU;
T_flow = ( U*UT + T*E*TU ) / ( U*UE + E*EU );

% simulated mean wage
Wmn = mean( wage );

% mean log wage
logwmn = mean( log(wage) );

% ub to wage
b_Wmn = b/Wmn;

% firing costs to wage
F_Wmn = F/Wmn;

% simulated mean output and productivity
Ymn = mean( ysim );
Ytot = Ymn*E;

% Separation tenure ratio (for calibration)

% yearly tenure
ten_yr = floor(ten/12);

if F==0

    % max yearly tenure level for ref group
    L = 15;

    % initialize
    sr_tn = zeros(L+1,1);

    for itn = 1:L+1
        xx = xsim(ten_yr==itn);
        n = histcounts(xx, xbins)';
        h = n / sum(n);
        sr_tn(itn) = h'*sp;
    end
end

% ratio
sr_tn = sr_tn(1) / mean(sr_tn(2:L+1));

% tenure adjusted perm-temp differential
% y = log(wage);
% L = length(wage);
% X = double([ones(L,1), perm, ten_yr]);
% coefs = (X'*X)^(-1)*(X'*y);
% Wdiff = coefs(2);

% table
agg_stat = table( U, T, UE, EU, b, sr_tn, wr, UP, PU, TP, TU, T_inflow, T_outflow, T_flow, Ytot, ...
    Ymn, Wmn, logwmn, b_Wmn, F_Wmn, theta, kkappa, xp, theta, wnew_mn );

if F>0
    agg_stat.xt = xt;
    agg_stat.xhat = xhat;
    agg_stat.xtilde = xtilde;
end

% welfare
if compute_welfare

    % welfare
    agg_stat.cV = mean(simulations_outcomes.cV);
    agg_stat.cJ = mean(simulations_outcomes.cJ);
    agg_stat.cU = wr;
    agg_stat.cTot = U*agg_stat.cU + (1-U)*(agg_stat.cV+agg_stat.cJ);

end

time.statistics = toc;
time.total = toc(timer);
eql.time = time;
equilibrium_outcomes = eql;

end

