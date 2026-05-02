% Script chính chạy mô phỏng Poisson-Schrödinger tự nhất[cite: 1]
clear; clc;

% --- Cấu hình thông số ---
L = 15; N = 500; % 15 nm[cite: 1]
z = linspace(0, L, N);
dz_nm = z(2) - z(1);
m_star = 0.067 * 9.109e-31; % GaAs[cite: 1]
T = 300; % 300 K[cite: 1]
Nd_2D = 1e17; % n-doped GaAs layer[cite: 1]

% Mật độ tạp chất 3D silicon tại trung tâm (z=0 trong hình 1)[cite: 1]
Nd_3D = zeros(1, N);
Nd_3D(round(N/2-5):round(N/2+5)) = Nd_2D / 1e-9; % Dày ~1nm

% Tiềm năng giam giữ octic anharmonic (Eq. 2)[cite: 1]
V0 = 228; k = 5; b1 = -2; b2 = 0.3;
z_shift = (z - L/2); 
V_conf = V0 * (b1*(z_shift/k).^2 + b2*(z_shift/k).^8);

% --- Vòng lặp hội tụ ---
VH = zeros(1, N); % Khởi tạo ban đầu[cite: 1]
EF_old = 0;
mixing = 0.2; % Để ổn định vòng lặp[cite: 1]

for iter = 1:100
    % 1. Giải Schrödinger[cite: 1]
    [Ei, psi] = solve_schrodinger(z, V_conf + VH, m_star);
    
    % 2. Tính n(z) và mức Fermi EF[cite: 1]
    [nz, EF_new] = calculate_density(Ei(1:15), psi(:,1:15), Nd_2D, dz_nm*1e-9, m_star, T);
    
    % 3. Cập nhật VH từ Poisson[cite: 1]
    VH_new = calculate_VH(z, nz, Nd_3D);
    
    % 4. Kiểm tra hội tụ (Eq. image_bfc945.png)[cite: 1]
    if abs(EF_new - EF_old) < 1e-3
        fprintf('Hội tụ thành công tại vòng lặp %d. EF = %.4f meV\n', iter, EF_new);
        break;
    end
    
    % Trộn thế năng để lặp bước tiếp theo[cite: 1]
    VH = (1-mixing)*VH + mixing*VH_new;
    EF_old = EF_new;
end

% --- Sau khi hội tụ, có thể gọi lệnh vẽ hình tại đây ---
% --- Lệnh vẽ hình cho Hình 2a (Nd = B = F = 0) ---

% 1. Thiết lập cửa sổ đồ thị
figure('Color', 'w', 'Name', 'Figure 2a: Nd=B=F=0');
hold on; grid on;

% 2. Vẽ tiềm năng giam giữ (Confining Potential)
% Chuyển z về hệ tọa độ nm để khớp với bài báo
z_plot = z - (L/2); 
plot(z_plot, V_conf, 'b', 'LineWidth', 2.5, 'DisplayName', 'Confining potential');

% 3. Vẽ 4 mức năng lượng đầu tiên và hàm mật độ xác suất tương ứng
colors = {'r', 'g', 'k', 'm'};
labels = {'Ground state', '1st excited', '2nd excited', '3rd excited'};
scale_factor = 500; % Hệ số để hiển thị hàm sóng trên thang năng lượng

for i = 1:4
    % Vẽ đường mức năng lượng Ei (đường thẳng đứt nét)[cite: 1]
    line([min(z_plot), max(z_plot)], [Ei(i), Ei(i)], ...
        'Color', colors{i}, 'LineStyle', '--', 'LineWidth', 1, 'HandleVisibility', 'off');
    
    % Tính và vẽ mật độ xác suất |psi|^2[cite: 1]
    % Hàm sóng được dịch lên mức năng lượng tương ứng để giống Hình 2a[cite: 1]
    prob_density = (psi(:,i).^2 / (dz_nm*1e-9)) * (dz_nm*1e-9 * scale_factor) + Ei(i);
    plot(z_plot, prob_density, 'Color', colors{i}, 'LineWidth', 2, 'DisplayName', labels{i});
end

% 4. Định dạng đồ thị theo tiêu chuẩn của bài báo[cite: 1]
xlabel('z (nm)', 'FontSize', 12, 'FontWeight', 'bold');
ylabel('Energy (meV)', 'FontSize', 12, 'FontWeight', 'bold');
title('Anharmonic GaAs QW (N_d = B = F = 0)', 'FontSize', 14);

% Giới hạn trục tọa độ tương ứng với Hình 2a trong tài liệu[cite: 1]
axis([-7.5 7.5 -500 1000]); 
legend('Location', 'northeast', 'FontSize', 10);

hold off;
