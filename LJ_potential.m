function V = LJ_potential(rs, rmin, ep, r_c)
    V_c = 4*ep*((2^(-1/6)*rmin/r_c)^12-(2^(-1/6)*rmin/r_c)^6);
    V = 4*ep*((2^(-1/6)*rmin./rs).^12-(2^(-1/6)*rmin./rs).^6) - V_c;
    V(rs > r_c) = 0;
end