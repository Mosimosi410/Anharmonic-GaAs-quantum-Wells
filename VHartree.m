function VH = calculate_VH(z, nz, Nd_3D)
    %nz: mật độ electron (m^-3) từ Eq. 9
    %Nd_3D: mật độ tạp chất silicon (m^-3)
    
    N = length(z);
    dz = (z(2) - z(1)) * 1e-9; % Chuyển nm sang m
    epsilon = 13.1 * 8.854e-12; % Hằng số điện môi GaAs
    e = 1.602e-19; 

    % Mật độ điện tích tổng cộng rho(z)
    rho = e * (Nd_3D - nz); % Eq. 8

    % Thiết lập ma trận Poisson (sai phân bậc hai)
    A = (diag(-2*ones(N,1)) + diag(ones(N-1,1),1) + diag(ones(N-1,1),-1)) / dz^2;
    
    % Điều kiện biên Dirichlet: VH = 0 tại hai đầu[cite: 1]
    A(1,:) = 0; A(1,1) = 1;
    A(N,:) = 0; A(N,N) = 1;

    % Vế phải phương trình Poisson: -rho/epsilon
    B = -rho / epsilon;
    B(1) = 0; B(N) = 0;

    % Giải hệ phương trình và chuyển sang meV
    VH_volts = A \ B';
    VH = VH_volts' * 1000; % Đơn vị meV[cite: 1]
end
