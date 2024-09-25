% M = 16;
% d = [0:M-1];
% y = qammod(d,M,'PlotConstellation',true);
% 
% % 設置坐標軸的字體大小
% ax = gca;
% ax.FontSize = 14;
% 
% M = 16;
% d = [0:M-1];
% y = qammod(d,M,'PlotConstellation',false);
% y_normalized = y / max(abs(y));
% catter(real(y_normalized), imag(y_normalized));

% M = 16;
% d = [0:M-1];
% y = qammod(d, M);
% 
% % 標準化星座點
% y_normalized = y / sqrt(mean(abs(y).^2));
% 
% % 繪製標準化後的星座圖
% scatterplot(y_normalized);
% scatter(real(y_normalized), imag(y_normalized));
% 
% grid on;

% 
% M = 16;
% d = [0:M-1];
% y = qammod(d, M);
% 
% % 標準化星座點
% y_normalized = y / sqrt(mean(abs(y).^2));
% 
% % 繪製標準化後的星座圖
% figure;
% scatter(real(y_normalized), imag(y_normalized), 'filled');
% axis square;
% grid on;
% 
% % 設置背景為白色
% set(gcf, 'Color', 'white');
% 
% % 設置縱橫軸的顏色
% ax = gca;
% ax.XColor = 'black';
% ax.YColor = 'black';

M = 16;
d = [0:M-1];
y = qammod(d, M, 'PlotConstellation', false);

% 標準化數據
y_normalized = y / sqrt(mean(abs(y).^2));

% 繪製標準化後的星座圖，背景為白色
figure;
scatterplot(y_normalized);
set(gcf, 'Color', 'white');

% 設置坐標軸的字體大小
ax = gca;
ax.FontSize = 14;

