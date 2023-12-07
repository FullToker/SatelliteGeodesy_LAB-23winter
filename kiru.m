% 读取数据文件
data_file_path = '/Users/lixipeng/Desktop/lab1_data/KIRU_ig1.xyz';
ref_file_path = '/Users/lixipeng/Desktop/lab1_data/ITRF2008_GNSS.SSC.txt';

% 使用importdata函数读取数据文件
data = importdata(data_file_path);

% 提取 MJD 和坐标
MJD = data.data(:, 2); % MJD，第二列
X = data.data(:, 5); % X坐标，第五列
Y = data.data(:, 6); % Y坐标，第六列
Z = data.data(:, 7); % Z坐标，第七列

% 我提取了我们需要的三个站点的坐标，不用读取，直接type in
ref_coords = [2251420.790, 862817.219, 5885476.705]; % KIRU站的坐标

% 坐标转换
local_coords = transform_to_local_horizontal(X, Y, Z, ref_coords);

% 日期转换公式
JD = MJD + 2400000.5; % 修正为儒略日
T = (JD - 2451545.0) / 36525; % 单位时间：世纪
Year = 2000 + T * 100; % 年份

% 最小二乘法拟合，向北（N）
omega = 2 * pi; % 一年一个周期
A = [ones(size(Year)), Year, cos(omega * Year), sin(omega * Year)];
lambda = 0; % 岭回归参数
beta_north = (A' * A + lambda * eye(size(A, 2))) \ (A' * local_coords(:, 3)); % 拟合北向（N）

% 最小二乘法拟合，向东（E）
beta_east = (A' * A + lambda * eye(size(A, 2))) \ (A' * local_coords(:, 2)); % 拟合东向（E）

% 最小二乘法拟合，高度（H）
beta_height = (A' * A + lambda * eye(size(A, 2))) \ (A' * local_coords(:, 1)); % 拟合高度（H）

% 计算线性趋势和季节性信号
linear_trend_north = A * beta_north;
linear_trend_east = A * beta_east;
linear_trend_height = A * beta_height;

% 计算残差
residuals_north = local_coords(:, 3) - linear_trend_north;
residuals_east = local_coords(:, 2) - linear_trend_east;
residuals_height = local_coords(:, 1) - linear_trend_height;

% 计算变化率（米/年）
rate_north = beta_north(2); % 北向变化率
rate_east = beta_east(2); % 东向变化率
rate_height = beta_height(2); % 高度变化率

% 北向（N）年周期信号的振幅
amplitude_north = sqrt(beta_north(3)^2 + beta_north(4)^2);

% 东向（E）年周期信号的振幅
amplitude_east = sqrt(beta_east(3)^2 + beta_east(4)^2);

% 高度（H）年周期信号的振幅
amplitude_height = sqrt(beta_height(3)^2 + beta_height(4)^2);

% 绘制原始数据、线性趋势和残差，向北（N）
figure;
subplot(3, 1, 1);
plot(Year, local_coords(:, 2), 'b');
title('原始北向(N)坐标数据');

subplot(3, 1, 2);
plot(Year, linear_trend_north, 'r');
title('北向(N)坐标的线性趋势');

subplot(3, 1, 3);
plot(Year, residuals_north, 'k');
title('北向(N)坐标的残差');
xlabel('时间(年)');
ylabel('坐标值(m)');

% 绘制原始数据、线性趋势和残差，向东（E）
figure;
subplot(3, 1, 1);
plot(Year, local_coords(:, 1), 'b');
title('原始东向(E)坐标数据');

subplot(3, 1, 2);
plot(Year, linear_trend_east, 'r');
title('东向(E)坐标的线性趋势');

subplot(3, 1, 3);
plot(Year, residuals_east, 'k');
title('东向(E)坐标的残差');
xlabel('时间(年)');
ylabel('坐标值(m)');

% 绘制原始数据、线性趋势和残差，高度（H）
figure;
subplot(3, 1, 1);
plot(Year, local_coords(:, 3), 'b');
title('原始高度(H)数据');

subplot(3, 1, 2);
plot(Year, linear_trend_height, 'r');
title('高度(H)的线性趋势');

subplot(3, 1, 3);
plot(Year, residuals_height, 'k');
title('高度(H)的残差');
xlabel('时间(年)');
ylabel('坐标值(m)');

% 绘制变化率，向北（N）
figure;
plot(Year, rate_north * Year + beta_north(1), 'b');
title('北向(N)坐标的变化率');
xlabel('时间(年)');
ylabel('变化率 (米/年)');

% 绘制变化率，向东（E）
figure;
plot(Year, rate_east * Year + beta_east(1), 'b');
title('东向(E)坐标的变化率');
xlabel('时间(年)');
ylabel('变化率 (米/年)');

% 绘制变化率，高度（H）
figure;
plot(Year, rate_height * Year + beta_height(1), 'b');
title('高度(H)的变化率');
xlabel('时间(年)');
ylabel('变化率 (米/年)');

function local_coords = transform_to_local_horizontal(X, Y, Z, ref_coords)
    % 计算参考点的地理纬度和经度
phi = atan2(ref_coords(2), ref_coords(1)); % 经度
lambda = atan2(ref_coords(3), sqrt(ref_coords(1)^2 + ref_coords(2)^2)); % 纬度

% 旋转矩阵
R2 = [cos(-lambda) 0 -sin(-lambda);
      0 1 0;
      sin(-lambda) 0 cos(-lambda)];

R3 = [cos(phi) sin(phi) 0;
      -sin(phi) cos(phi) 0;
      0 0 1];

    % 从地心地固坐标到局部坐标的转换
    local_coords = zeros(length(X), 3); % 初始化局部坐标数组
    for i = 1:length(X)
        % 原始坐标减去参考点坐标
        coord = [X(i) - ref_coords(1); Y(i) - ref_coords(2); Z(i) - ref_coords(3)];
        % 应用旋转矩阵
        local_coord = R2 * R3 * coord;
        local_coords(i, :) = local_coord';
    end
end
