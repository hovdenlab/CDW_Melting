%% Setup MC Simulation
mul = 3; % particle count parameter
meltIter = 250; % initial melting iterations
relaxIter = 50; % iterations for each site fluctuation
temp_idxs = 0:0.25:15; % arb. temp slices
num_add = 50; % maximum number of possible sites added at each temperature
num_remove = 150; % maximum number of possible sites removed at each temperature

% LJ parameters
rmin = 1; % location of well minimum
r_c = 15*rmin; % potential cutoff
ep = 10; % potential well depth (arb.)

%% Tools
getR = @(x) sqrt( x(:,1).^2 + x(:,2).^2 );

rng('shuffle')
sd = rng; 

%% Grid Setup
numGridLT = 800;   m = (-numGridLT:numGridLT);
[m,n] = meshgrid(m,m);
numGridSL = floor(numGridLT/2);  j = (-numGridSL:numGridSL);
[j,k] = meshgrid(j,j);

m = m(:); n = n(:); j = j(:); k = k(:);

th = 0;
a = 1 * 1/sqrt(13);
a1 = a*[cosd(th), sind(th)]; a2 = a*[cosd(th+120),sind(th+120)];
aSL = sqrt(13)*a;
SL1 = aSL*[cosd(0), sind(0)]; SL2 = aSL*[cosd(120), sind(120)];

sg = aSL/2;

posLT = m*a1 + n*a2;
posSL = j*SL1 + k*SL2;

xCutoffSL = aSL*21; yCutoffSL = xCutoffSL/7*4*sqrt(3);
xCutoffSL = mul*xCutoffSL; yCutoffSL = mul*yCutoffSL;

rd = 4;
a1 = round(a1,rd); a2 = round(a2,rd);
aSL = round(aSL,rd); SL1 = round(SL1, rd); SL2 = round(SL2, rd);
posLT = round( posLT,rd); posSL = round( posSL, rd);
xCutoffSL = round( xCutoffSL, rd); yCutoffSL = round( yCutoffSL, rd);

posSL( posSL(:,1) >=  xCutoffSL ,: ) = [];
posSL( posSL(:,1) <   0         ,: ) = [];
posSL( posSL(:,2) >=  yCutoffSL ,: ) = [];
posSL( posSL(:,2) <   0         ,: ) = [];

rposLT = sqrt( (posLT(:,1)-xCutoffSL/2).^2 + (posLT(:,2)-yCutoffSL/2).^2);

posLT( rposLT > xCutoffSL/2*0.7, : ) = [];

numSiteLT = length(posLT);
numSiteSL = length(posSL);


%% Monte Carlo Setup
% Iteration Setup
kBTs = temp_idxs;
numT = length(kBTs);
iter = meltIter*ones([numT,1]);
iter(1) = 0;

% Number of sites to remove
iter_remove = num_remove*ones([numT,1]);
iter_remove(1) = 0;

% Number of sites to add
iter_add = num_remove*ones([numT,1]);
iter_add(1) = 0;

% Random Deviation Limit
delta = aSL*0.1;

posLT_PLD = zeros([size(posLT),numT]);
posSLs    = -1*ones([length(posSL)+num_add,2,numT]);
pldScal = 0.1;

E_iters = zeros(meltIter,numT);
E_fluc_iters = zeros(num_remove+1+num_add, numT);

numSiteSLs = numSiteSL*ones(numT, 1);
delta_sites = zeros(numT, 1);

%% Iterate
for indT = 1:numT
    kBT = kBTs(indT);
    disp(kBT)

    if indT == 1 % skip temp = 0
        posSLs(1:length(posSL),:,indT)    = posSL;
        posLT_PLD(:,:,indT) = calcPLD_gaussian(posLT, posSL, aSL/2, aSL/2, pldScal);
        continue
    end

    % Initial melting of charge lattice
    for ct = 1:iter(indT)
        for indSL = randperm(numSiteSL)
            pos_cur = posSL(indSL,:);
            dPos = round(delta*(2*rand(1,2)-1),rd);
            
            dR_nns = pos_cur - posSL;
            dR_nns(:,1) = dR_nns(:,1) - round(dR_nns(:,1) / xCutoffSL) * xCutoffSL;
            dR_nns(:,2) = dR_nns(:,2) - round(dR_nns(:,2) / yCutoffSL) * yCutoffSL;
            
            r_nns = getR(dR_nns);
            dR_nns( r_nns > r_c | r_nns <= 0, :) = [];
            dR_nns_new = dR_nns + dPos;
            r_nns_cur = getR(dR_nns);
            r_nns_new = getR(dR_nns_new);

            curE = sum(LJ_potential(r_nns_cur, rmin, ep, r_c));
            newE = sum(LJ_potential(r_nns_new, rmin, ep, r_c));
            dE = newE - curE;

            E_iters(ct, indT) = E_iters(ct, indT) + curE;

            if rand() < exp(-dE/(kBT))
                pos_new = pos_cur + dPos;
                posSL(indSL,1) = pos_new(1) - floor(pos_new(1)/yCutoffSL)* xCutoffSL;
                posSL(indSL,2) = pos_new(2) - floor(pos_new(2)/xCutoffSL)* yCutoffSL;
                E_iters(ct, indT) = E_iters(ct, indT) + dE;
            end
        end
    end

    [posSL, E_ix] = melt_all_iter(posSL, kBT, rmin, r_c, ep, relaxIter, xCutoffSL, yCutoffSL);
    E_fluc_iters(num_remove+1, indT) = E_ix(end);

    % Sweep add/remove sites and find minimum energy state
    % Remove up to num_remove sites, and iterate relaxIter times
    kill_posSL = posSL;
    kill_posSLs = -1*ones([size(kill_posSL), iter_remove(indT)]);
    for ix = iter_remove(indT):-1:1

        kill_options = randsample(length(kill_posSL), 100);
        psi6_options = zeros(size(kill_options));
        for ko_ix = 1:length(kill_options)
            psi6_options(ko_ix) = psi6_6nn(kill_posSL(kill_options(ko_ix),:), kill_posSL);
        end
        [mpsi6, mp_ix] = mink(psi6_options, 1);
        kill_ix = kill_options(mp_ix);
    
        kx = kill_posSL(kill_ix, 1);
        ky = kill_posSL(kill_ix, 2);
        kill_posSL(kill_ix,:) = [];
        
        [kill_posSL, E_ix] = melt_all_iter(kill_posSL, kBT, rmin, r_c, ep, relaxIter, xCutoffSL, yCutoffSL);
        E_fluc_iters(ix, indT) = E_ix(end);

        kill_posSLs(1:length(kill_posSL),:,ix) = kill_posSL;

    end

    % Add up to num_add sites, and iterate relaxIter times
    add_posSL = posSL;
    add_posSLs = -1*ones([length(add_posSL)+num_add, 2, iter_remove(indT)]);
    for ix = 1:iter_add(indT)

        xs = add_posSL(:,1);
        ys = add_posSL(:,2);

        x_options = datasample(2:xCutoffSL-2,100);
        y_options = datasample(2:yCutoffSL-2,100);
        
        density_options = zeros(size(x_options));
        for add_ix = 1:length(x_options)
            numcount = length(find(xs>=(x_options(add_ix)-1) & xs<=(x_options(add_ix)+1) & ys>=(y_options(add_ix)-1) & ys<=(y_options(add_ix)+1)));
            density_options(add_ix) = numcount;
        end
        [~, mp_ix] = mink(density_options, 1);
        add_x = x_options(mp_ix);
        add_y = y_options(mp_ix);
    
        add_posSL = [add_posSL; [add_x add_y]];
        
        [add_posSL, E_ix] = melt_all_iter(add_posSL, kBT, rmin, r_c, ep, relaxIter, xCutoffSL, yCutoffSL);
        E_fluc_iters(num_remove+1+ix, indT) = E_ix(end);

        add_posSLs(1:length(add_posSL),:,ix) = add_posSL;

    end

    % Select lowest energy configuration
    options = 1:length(E_fluc_iters(:,indT));
    Pfit = polyfit(options, E_fluc_iters(:,indT), 2);
    m_ix = round(-Pfit(2)/(2*Pfit(1)));
    if m_ix < 1
        m_ix = 1;
    end
    if m_ix > num_remove+1+num_add
        m_ix = num_remove+1+num_add;
    end

    delta_sites(indT) = m_ix - (num_remove+1);
    numSiteSL = numSiteSL + delta_sites(indT);
    numSiteSLs(indT) = numSiteSL;

    if m_ix <= num_remove
        posSL = kill_posSLs(1:numSiteSL,:,m_ix);
    elseif m_ix >= num_remove+2
        posSL = add_posSLs(1:numSiteSL,:,m_ix-(num_remove+1));
    end
 
    posSLs(1:length(posSL),:,indT)    = posSL;
    posLT_PLD(:,:,indT) = calcPLD_gaussian(posLT, posSL, aSL/2, aSL/2, pldScal);
    
end

save( sprintf('data/melt_fluctuate_charge_%d.mat',(sd.Seed)) )
return
