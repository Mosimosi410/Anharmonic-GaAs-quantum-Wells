% Chương trình vẽ Hình 2a: Trạng thái Nd = B = F = 0
% Giải phương trình Schrödinger cho giếng lượng tử octic anharmonic

clear; clc;

%% 1. Khai báo hằng số vật lý và thông số cấu trúc
hbar = 1.05457e-34;         % Hằng số Planck rút gọn (J.s)
e = 1.60218e-19;            % Điện tích electron (C)
m0 = 9.10938e-31;           % Khối lượng nghỉ electron (kg)
m_star = 0.067 * m0;        % Khối lượng hiệu dụng GaAs

L = 15;                     % Chiều dài cấu trúc (nm)
N = 1000;                   % Số điểm lưới để đạt độ chính xác cao
z = linspace(-L/2, L/2, N); % Tọa độ z từ -7.5 đến 7.5 nm
dz = (z(2) - z(1)) * 1e-9;  % Bước lưới (m)

%% 2. Thiết lập tiềm năng giam giữ V_conf (Octic Anharmonic)
V0 = 228;                   % Độ sâu giếng (meV)
k_p = 5;                    % Tham số k (nm)[cite: 1]
beta1 = -2;                 % Hệ số beta1[cite: 1]
beta2 = 0.3;                % Hệ số beta2[cite: 1]

% Công thức: V_conf(z) = V0 * [beta1*(z/k)^2 + beta2*(z/k)^8][cite: 1]
V_conf = V0 * (beta1 * (z/k_p).^2 + beta2 * (z/k_p).^8); 

%% 3. Xây dựng ma trận Hamiltonian (H) bằng phương pháp sai phân hữu hạn (FDM)
% Tính hệ số động năng t (đơn vị meV)[cite: 1]
t = (hbar^2 / (2 * m_star * dz^2)) / (e * 1e-3); 

% Tạo ma trận tridiagonal (Eq. 13)[cite: 1]
main_diag = 2*t + V_conf;           
off_diag  = -t * ones(1, N-1);      

H = diag(main_diag) + diag(off_diag, 1) + diag(off_diag, -1);

%% 4. Giải bài toán trị riêng để tìm mức năng lượng (Ei) và hàm sóng (psi)
[psi, E] = eig(H);
E_levels = diag(E);                 % Các giá trị năng lượng (meV)[cite: 1]

%% 5. Vẽ đồ thị mô phỏng Hình 2a
figure('Color', 'w', 'Name', 'Figure 2a: Nd=B=F=0');
hold on; grid on;

% Vẽ tiềm năng giam giữ (đường màu xanh lam đậm)[cite: 1]
plot(z, V_conf, 'b', 'LineWidth', 2.5, 'DisplayName', 'Confining potential');

% Thiết lập màu sắc và nhãn cho 4 mức năng lượng đầu tiên[cite: 1]
colors = {'r', 'g', 'k', 'm'};
labels = {'Ground state', '1st excited', '2nd excited', '3rd excited'};
scale = 500; % Hệ số tỉ lệ để hiển thị xác suất trên trục năng lượng

for i = 1:4
    % Vẽ đường mức năng lượng Ei (đường thẳng đứt nét)
    line([-L/2, L/2], [E_levels(i), E_levels(i)], 'Color', colors{i}, 'LineStyle', '--', 'LineWidth', 1);
    
    % Vẽ mật độ xác suất |psi|^2 dịch chuyển theo mức năng lượng tương ứng[cite: 1]
    % Chú ý: psi(:,i) cần được chuẩn hóa để diện tích dưới đường psi^2 bằng 1
    prob_density = (psi(:,i).^2 / dz) * (dz * scale) + E_levels(i);
    plot(z, prob_density, colors{i}, 'LineWidth', 1.8, 'DisplayName', labels{i});
end

% Định dạng trục và chú thích giống bài báo[cite: 1]
xlabel('z (nm)', 'FontSize', 12);
ylabel('Energy (meV)', 'FontSize', 12);
title('Anharmonic GaAs QW (N_d = B = F = 0)', 'FontSize', 14);
axis([-7.5 7.5 -500 1000]); % Giới hạn trục theo Hình 2a[cite: 1]
legend('Location', 'northeast');

% Hiển thị giá trị năng lượng cụ thể trên Command Window
fprintf('Các mức năng lượng tính toán được (meV):\n');
for i = 1:4
    fprintf('%s: %.2f meV\n', labels{i}, E_levels(i));
end
