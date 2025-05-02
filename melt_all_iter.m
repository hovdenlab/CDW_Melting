function [posSL, E_iters] = melt_all_iter(posSL, temp, rmin, r_c, ep, numIter, xCutoffSL, yCutoffSL)    
    getR = @(x) sqrt( x(:,1).^2 + x(:,2).^2 );
    rd = 4;
   
    numSiteSL = length(posSL);
    
    % Random Deviation Limit
    aSL = 1;
    delta = aSL*0.1;
    
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
            E_iters(ct) = E_iters(ct) + curE;

            dE = newE - curE;

            if rand() < exp(-dE/temp)
                pos_new = pos_cur + dPos;
                posSL(indSL,1) = pos_new(1) - floor(pos_new(1)/yCutoffSL)* xCutoffSL;
                posSL(indSL,2) = pos_new(2) - floor(pos_new(2)/xCutoffSL)* yCutoffSL;
                E_iters(ct) = E_iters(ct) + dE;
            end
        end 
    end
end

