function [posSL, scalings, E_iters, xCutoffSL, yCutoffSL] = melt_single_iter(posSL, temp, rmin, r_c, ep, numIter, xCutoffSL, yCutoffSL)
    getR = @(x) sqrt( x(:,1).^2 + x(:,2).^2 );
    rd = 4;
   
    numSiteSL = length(posSL);
    
    % Random Deviation Limit
    aSL = 1;
    delta = aSL*0.1;
    dScale = 0.01;
    
    scalings = ones(numIter, 1);
    E_iters = zeros(numIter, 1);
    for ct = 1:numIter       
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

            newE = sum(LJ_potential(r_nns_new, rmin, ep, r_c));
            curE = sum(LJ_potential(r_nns_cur, rmin, ep, r_c));
            dE = newE - curE;

            if rand() < exp(-dE/temp)
                pos_new = pos_cur + dPos;
                posSL(indSL,1) = pos_new(1) - floor(pos_new(1)/yCutoffSL)* xCutoffSL;
                posSL(indSL,2) = pos_new(2) - floor(pos_new(2)/xCutoffSL)* yCutoffSL;
            end
        end 

        xs = posSL(:,1);
        ys = posSL(:,2);
        rs = allRsPBC(xs, ys, r_c, xCutoffSL, yCutoffSL);
        curE = sum(LJ_potential(rs, rmin, ep, r_c));
        E_iters(ct) = curE;

        scale = 1+dScale*(2*rand()-1);
        rs = allRsPBC(xs*scale, ys*scale, r_c, xCutoffSL*scale, yCutoffSL*scale);
        newE = sum(LJ_potential(rs, rmin, ep, r_c));

        dE = newE-curE;
        if rand() < exp(-dE/(temp))
            E_iters(ct) = newE;
            scalings(ct) = scale;
            posSL = posSL*scale;
            xCutoffSL = xCutoffSL*scale;
            yCutoffSL = yCutoffSL*scale;
        end

    end
end

