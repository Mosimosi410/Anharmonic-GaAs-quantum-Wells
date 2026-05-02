function [nz, EF] = calculate_density(Ei, psi, Nd_2D, dz_m, m_star, T)
    % Ei: các mức năng lượng (meV)[cite: 1]
    % psi: hàm sóng đã giải từ Schrodinger[cite: 1]
    
    kB = 1.38e-23;
    hbar = 1.054e-34;
    e = 1.602e-19;
    
    % Hệ số mật độ trạng thái 2D: (m* kB T) / (pi hbar^2)[cite: 1]
    pref = (m_star * kB * T) / (pi * hbar^2);
    
    % Tìm mức Fermi EF sao cho thỏa mãn neutrality condition (Eq. 10)[cite: 1]
    % f(EF) = Nd_2D - sum( Ni(EF) ) = 0
    f_neutrality = @(EF_val) Nd_2D - sum(pref * log(1 + exp((EF_val - Ei) * e * 1e-3 / (kB * T))));
    
    % Tìm nghiệm EF (bắt đầu từ mức năng lượng thấp nhất)
    EF = fzero(f_neutrality, Ei(1)); 

    % Tính mật độ electron n(z) theo Eq. 9[cite: 1]
    N = size(psi, 1);
    nz = zeros(1, N);
    for i = 1:size(psi, 2) % Duyệt qua các mức năng lượng
        occupancy = pref * log(1 + exp((EF - Ei(i)) * e * 1e-3 / (kB * T)));
        nz = nz + occupancy * (psi(:,i).^2 / dz_m)'; 
    end
end
