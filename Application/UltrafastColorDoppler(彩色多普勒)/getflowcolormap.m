%% 自定义colormap

function customMap = getflowcolormap()
numColors = 256;

selectedColors =[
    44,234,255;
    6,208,255;
    0,171,253;
    0,127,254;
    1,88,255;
    0,47,249;
    0,18,216;
    0,5,123;
    29, 0, 0;
    161,4,11;
    203,17,22;
    227,43,31;
    244,81,40;
    255,127,51;
    255,170,61;
    255,208,70;
    244,230,133
    ];

numPoints = size(selectedColors,1);

xq = linspace(1, numPoints, numColors);
r = interp1(1:numPoints, selectedColors(:, 1), xq, 'pchip');
g = interp1(1:numPoints, selectedColors(:, 2), xq, 'pchip');
b = interp1(1:numPoints, selectedColors(:, 3), xq, 'pchip');


% 创建colormap矩阵
customMap = [r(:), g(:), b(:)]/256;


% data = rand(100);
% imagesc(data); colormap(customMap); caxis([0 1]); colorbar;
% 
% 
% [X, Y, Z] = peaks(50);
% surf(X, Y, Z, 'EdgeColor', 'none');
% colormap(customMap);
% caxis([-10 10]);
% colorbar;

end

