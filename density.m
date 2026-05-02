clear; clc;

%% 1. Thông số vật lý và thiết lập lưới
m_star = 0.067 * 9.109e-31; % Khối lượng hiệu dụng GaAs
hbar = 1.054e-34;
e = 1.602e-19;
kB = 1.38e-23;
T = 300;                    % Nhiệt độ (K)[cite: 1]
epsilon = 13.1 * 8.854e-12; % hằng số điện môi GaAs[cite: 1]

L = 15e-9;                  % Chiều dài 15 nm[cite: 1]
N = 500;                    % Số điểm lưới
z = linspace(0, L, N);
dz = z(2) - z(1);

%% 2. Thiết lập nồng độ pha tạp Nd(z)[cite: 1]
% Lớp silicon pha tạp dày 2nm ở giữa (z=7.5nm)[cite: 1]
Nd_2D = 1e17;               % Nồng độ 2D (m^-2)[cite: 1]
Nd_3D = zeros(1, N);
center_idx = round(N/2);
width_idx = round(1e-9 / dz); % 1nm mỗi bên trung tâm
Nd_3D(center_idx-width_idx : center_idx+width_idx) = Nd_2D / 2e-9; 

%% 3. Vòng lặp giải tự nhất (Self-consistent Loop)[cite: 1]
VH = zeros(1, N);           % Khởi tạo thế Hartree ban đầu[cite: 1]
EF_old = 0;
tolerance = 1e-3;           % Sai số 10^-3 meV[cite: 1]
max_iter = 100;
converged = false;

% Tính thế năng giam giữ gốc V_conf[cite: 1]
V0 = 228; k = 5; b1 = -2; b2 = 0.3;
z_nm = (z - L/2) * 1e9; % dịch tâm về 0 cho công thức potential
V_conf = V0 * (b1*(z_nm/k).^2 + b2*(z_nm/k).^8); 

for iter = 1:max_iter
    % --- BƯỚC 1: GIẢI PHƯƠNG TRÌNH SCHRÖDINGER (Eq. 11 & 13)[cite: 1] ---
    V_total = V_conf + VH; 
    t = (hbar^2 / (2 * m_star * dz^2)) / (e * 1e-3); % meV
    H = diag(2*t + V_total) + diag(-t*ones(1,N-1), 1) + diag(-t*ones(1,N-1), -1);
    [psi, E] = eig(H);
    Ei = diag(E); 
    
    % --- BƯỚC 2: TÌM MỨC FERMI EF (Neutrality Condition - Eq. 10)[cite: 1] ---
    % Giải phương trình: Nd = sum( (m*kB*T)/(pi*hbar^2) * log(1 + exp((EF-Ei)/kBT)) )
    const_pref = (m_star * kB * T) / (pi * hbar^2);
    f_EF = @(EF) Nd_2D - sum(const_pref * log(1 + exp((EF - Ei) * e * 1e-3 / (kB * T))));
    EF_new = fzero(f_EF, Ei(1)); % Tìm EF sao cho trung hòa điện[cite: 1]
    
    % --- BƯỚC 3: TÍNH MẬT ĐỘ ELECTRON n(z) (Eq. 9)[cite: 1] ---
    nz = zeros(1, N);
    for i = 1:10 % Xét 10 mức năng lượng đầu tiên
        occupancy = const_pref * log(1 + exp((EF_new - Ei(i)) * e * 1e-3 / (kB * T)));
        nz = nz + occupancy * (psi(:,i).^2 / dz)'; 
    end
    
    % --- BƯỚC 4: GIẢI PHƯƠNG TRÌNH POISSON TÌM VH(new) (Eq. 8)[cite: 1] ---
    rho = e * (Nd_3D - nz);
    A_poisson = (diag(-2*ones(1,N)) + diag(ones(1,N-1),1) + diag(ones(1,N-1),-1)) / dz^2;
    A_poisson(1,:) = 0; A_poisson(1,1) = 1; % Biên VH=0[cite: 1]
    A_poisson(N,:) = 0; A_poisson(N,N) = 1;
    B_poisson = -rho / epsilon; B_poisson(1) = 0; B_poisson(N) = 0;
    VH_new = (A_poisson \ B_poisson')' * (1000/e); % meV
    
    % --- KIỂM TRA HỘI TỤ (Ef_new - Ef_old) < 10^-3 meV[cite: 1] ---
    if abs(EF_new - EF_old) < tolerance
        fprintf('Hội tụ tại vòng lặp % d. EF = %.4f meV\n', iter, EF_new);
        converged = true;
        break;
    end
    
    % Cập nhật (sử dụng mixing để ổn định vòng lặp)
    mixing = 0.2;
    VH = (1 - mixing) * VH + mixing * VH_new;
    EF_old = EF_new;
end
