function [rs, thetas] = rtheta(xs, ys, rmax)
    num = length(xs);
    rs = [];
    thetas = [];
    for i = 1:num
        for j = i+1:num
            r =  sqrt((xs(i)-xs(j))^2 + (ys(i)-ys(j))^2);
            if (r <= rmax)
                rs = [rs r];
                theta = atan2((ys(i)-ys(j)),(xs(i)-xs(j)));
                thetas = [thetas theta];
            end
        end
    end
end