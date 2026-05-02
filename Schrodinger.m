function [Ei, psi] = solve_schrodinger(z, V_total, m_star)
    % V_total: tổng thế năng (meV) = V_conf + VH[cite: 1]
    
    N = length(z);
    dz = (z(2) - z(1)) * 1e-9;
    hbar = 1.054e-34;
    e = 1.602e-19;
    
    % Hệ số động năng t (meV)[cite: 1]
    t = (hbar^2 / (2 * m_star * dz^2)) / (e * 1e-3);
    
    % Xây dựng ma trận Hamiltonian (Eq. 13)[cite: 1]
    main_diag = 2*t + V_total;
    off_diag = -t * ones(1, N-1);
    H = diag(main_diag) + diag(off_diag, 1) + diag(off_diag, -1);
    
    % Giải trị riêng
    [psi, E] = eig(H);
    Ei = diag(E); % Trả về vector các mức năng lượng (meV)[cite: 1]
end
