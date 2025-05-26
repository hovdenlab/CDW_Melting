fname = "test";
f = load("data/"+fname+".mat");

kill_iters = 100;
numIter = 150;

t_ix = 5;

idxs = f.posSLs;
temp_vals = f.temp_vals;
rmin = f.rmin;
r_c = f.r_c;
ep = f.ep;

posSL = idxs(:,:,t_ix); % get position at temp slice
xs = posSL(:,1);
ys = posSL(:,2);
temp = temp_vals(t_ix); % get temp
xCutoffSL0 = f.xCutoffSLs(1); % get original state
yCutoffSL0 = f.yCutoffSLs(1); % get original state

xCutoffSL = f.xCutoffSLs(t_ix); % get expanded state
yCutoffSL = f.yCutoffSLs(t_ix); % get expanded state

E_iters = zeros(numIter, kill_iters+1);
posSL_post = zeros(length(posSL), 2, kill_iters+1);

xCutoffSLs = zeros(kill_iters+1,1);
yCutoffSLs = zeros(kill_iters+1,1);

kill_pos = zeros(kill_iters, 2);
scalings = zeros(numIter, kill_iters);

% ensure energy convergence first
disp("Start Initial Energy Convergence")
[posSL, scaling_initial, E_initial_convergence, xCutoffSL, yCutoffSL] = melt_single_iter(posSL, temp, rmin, r_c, ep, numIter, xCutoffSL, yCutoffSL);

initial_pos = posSL;
posSL_post(:,:,1) = posSL;
E_iters(:,1) = E_initial_convergence;
xCutoffSLs(1) = xCutoffSL;
yCutoffSLs(1) = yCutoffSL;

nsample = 100;
for ix = 1:kill_iters
    disp("remove site #" + num2str(ix))
    
    % randomly select nsample sites and pick most disordered site
    kill_options = randsample(length(posSL), nsample);
    psi6_options = zeros(size(kill_options));
    for ko_ix = 1:length(kill_options)
        psi6_options(ko_ix) = psi6_6nn(posSL(kill_options(ko_ix),:), posSL);
    end
    [mpsi6, mp_ix] = min(psi6_options);
    kill_ix = kill_options(mp_ix);
    
    kx = posSL(kill_ix, 1);
    ky = posSL(kill_ix, 2);
    posSL(kill_ix,:) = [];

    % melt
    [posSL, scalings_iter, E_ix_iter, xCutoffSL, yCutoffSL] = melt_single_iter(posSL, temp, rmin, r_c, ep, numIter, xCutoffSL, yCutoffSL);

    E_iters(:,ix+1) = E_ix_iter;
    posSL_post(1:end-ix,:,ix+1) = posSL;
    kill_pos(ix,1) = kx;
    kill_pos(ix,2) = ky;
    scalings(:,ix) = scalings_iter;
    xCutoffSLs(ix+1) = xCutoffSL;
    yCutoffSLs(ix+1) = yCutoffSL;

end

[~, select_ix] = min(abs(xCutoffSLs - xCutoffSL0));
posSL = posSL_post(1:end-select_ix,:,select_ix);
num_sites_removed = select_ix;

save("data/" + fname + "_slice_" + num2str(t_ix) + ".mat", "temp_vals", "t_ix", "posSL", "num_sites_removed", "posSL_post", "kill_pos", "scalings", "xCutoffSLs", "yCutoffSLs", "E_iters", "numIter", "kill_iters")
