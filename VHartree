function VH = calculate_hartree(z, n_z, Nd_z)
    % Hàm tính thế năng Hartree bằng cách giải phương trình Poisson
    % z: vector vị trí (nm)
    % n_z: mật độ electron (m^-3)
    % Nd_z: mật độ tạp chất silicon (m^-3)

    N = length(z);
    dz = (z(2) - z(1)) * 1e-9; % Chuyển nm sang m
    
    % Các hằng số vật lý
    epsilon_0 = 8.854e-12;     % Hằng số điện môi chân không
    epsilon_r = 13.1;          % Hằng số điện môi của GaAs
    e = 1.602e-19;             % Điện tích nguyên tố
    
    % Thiết lập ma trận Poisson (đạo hàm bậc hai)
    % d^2(VH)/dz^2 = (e^2 / (eps0 * epsr)) * [Nd(z) - n(z)]
    
    % Ma trận sai phân bậc hai (Laplacian)
    A = (diag(-2 * ones(N, 1)) + diag(ones(N-1, 1), 1) + diag(ones(N-1, 1), -1)) / (dz^2);
    
    % Điều kiện biên: VH = 0 tại hai đầu cấu trúc (z=0 và z=L)
    A(1, :) = 0; A(1, 1) = 1;
    A(N, :) = 0; A(N, N) = 1;
    
    % Vế phải của phương trình (B)
    % Nhân thêm e để chuyển đơn vị sang eV hoặc meV tùy mục đích
    % Ở đây ta tính VH theo đơn vị meV để phù hợp với đồ thị bài báo
    B = (e / (epsilon_0 * epsilon_r)) * (Nd_z - n_z) * (1000); 
    
    % Áp dụng điều kiện biên cho vế phải
    B(1) = 0;
    B(N) = 0;
    
    % Giải hệ phương trình tuyến tính A * VH = B
    VH = A \ B'; 
    VH = VH'; % Trả về vector dòng
end
