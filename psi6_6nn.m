function psi6 = psi6_6nn(pos,neighbors)
        px = pos(1);
        py = pos(2);
        nxs = neighbors(:,1);
        nys = neighbors(:,2);
        dxs = nxs-px;
        dys = nys-py;
        [~, min_ix] = mink(dxs.^2+dys.^2, 7);
        theta = atan2(dys(min_ix(2:7)), dxs(min_ix(2:7)));
        psi6 = abs(mean(exp(6i*theta)));
end