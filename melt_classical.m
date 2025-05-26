%% Mode
mul = 2;
numIter = 10^3;
temp_vals = 0:1:10;

% LJ parameters
rmin = 1;
r_c = 15*rmin;
ep = 10;

save_mat_name = "test";

%% Tools
getR = @(x) sqrt( x(:,1).^2 + x(:,2).^2 );

boolDraw = false;

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
numT = length(temp_vals);
iter = numIter*ones([numT,1]);
iter(1) = 0;

% Random Deviation Limit
delta = aSL*0.1;

posLT_PLD = zeros([size(posLT),numT]);
posSLs    = zeros([size(posSL),numT]);
pldScal = 0.1;
E_iters = zeros(numIter,numT);

scale_freq = 25;
scalings = zeros(ceil(numIter/scale_freq),numT);

xCutoffSLs = zeros(numT,1);
yCutoffSLs = zeros(numT,1);

%% Relax Initial

side_scalings = 0.95:0.01:1.05;
E_scalings = zeros(size(side_scalings));

for s_ix = 1:length(side_scalings)    
    xs = posSL(:,1)*side_scalings(s_ix);
    ys = posSL(:,2)*side_scalings(s_ix);

    [rs, ~] = rtheta(xs, ys, r_c);
    E_scalings(s_ix) = sum(LJ_potential(rs, rmin, ep, r_c));
end

[min_E, m_ix] = min(E_scalings);
initial_scaling = side_scalings(m_ix);

disp("Initial Scaling " + num2str(initial_scaling))
     
posSL = posSL*initial_scaling;
xCutoffSL = xCutoffSL*initial_scaling;
yCutoffSL = yCutoffSL*initial_scaling;

xCutoffSLs(1) = xCutoffSL;
yCutoffSLs(1) = yCutoffSL;


%% Iterate
for indT = 1:numT
    kBT = temp_vals(indT);
    disp("temp: " + num2str(kBT))

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

            if rand() < exp(-dE/(kBT))
                pos_new = pos_cur + dPos;
                posSL(indSL,1) = pos_new(1) - floor(pos_new(1)/yCutoffSL)* xCutoffSL;
                posSL(indSL,2) = pos_new(2) - floor(pos_new(2)/xCutoffSL)* yCutoffSL;
            end
        end
        
        if mod(ct, scale_freq) == 0
            side_scalings = 0.9:0.01:1.1;
            E_scalings = zeros(size(side_scalings));
        
            for s_ix = 1:length(side_scalings)    
                xs = posSL(:,1)*side_scalings(s_ix);
                ys = posSL(:,2)*side_scalings(s_ix);
        
                [rs, ~] = rtheta(xs, ys, r_c);
                E_scalings(s_ix) = sum(LJ_potential(rs, rmin, ep, r_c));
            end
        
            [min_E, m_ix] = min(E_scalings);
            scalings(ct/scale_freq, indT) = side_scalings(m_ix);

            % disp("Scaling " + num2str(scalings(ct/scale_freq, indT)))
                 
            posSL = posSL*side_scalings(m_ix);
            xCutoffSL = xCutoffSL*side_scalings(m_ix);
            yCutoffSL = yCutoffSL*side_scalings(m_ix);
        end

        E_iters(ct, indT) = min_E;

    end
    posSLs(:,:,indT)    = posSL;
    posLT_PLD(:,:,indT) = calcPLD_gaussian(posLT, posSL, aSL/2, aSL/2, pldScal);
    xCutoffSLs(indT) = xCutoffSL;
    yCutoffSLs(indT) = yCutoffSL;
end

save("data/"+save_mat_name+".mat", "temp_vals", "posSLs", "posLT_PLD", "mul", "numIter", "rmin", "r_c", "ep", "scalings", "xCutoffSLs", "yCutoffSLs", "E_iters")

return

