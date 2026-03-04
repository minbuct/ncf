clc; clear; close all;
lamda = 0.3; 
stop_value = 1e-7;
ref_path = '数据集/椭球1000点(高噪声)/PointCloud_Standard.txt';
true_data = importdata(ref_path);
true_model = ellipsoidfit_leastsquares(true_data);
for i = 1:1:1
    niose_level = {'G0.35'};
    noise_name = char(niose_level(i));
    for j = 1:1:10
        instance_path = ['数据集/椭球1000点(高噪声)/', noise_name, 'num', num2str(j), '.txt'];
        data = importdata(instance_path); 
        model = Ours(data, lamda, stop_value); 
    end
end
figure('Color', 'w');
plot_ellipsoid(model, data); 
hold on;
scatter3(data(:,1), data(:,2), data(:,3), '.b'); 
axis equal; 
legend('Fitted Ellipsoid', 'Data points', 'Location', 'northwest');
set(gca, 'FontName', 'Times New Roman', 'FontSize', 14);
title(['Result of last case: ', noise_name, ' num', num2str(j)]);
