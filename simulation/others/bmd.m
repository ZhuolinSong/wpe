fit1 = cell2mat(struct2cell(load('bmd_fit1.mat')));
fit2 = cell2mat(struct2cell(load('bmd_fit2.mat')));
fit3 = cell2mat(struct2cell(load('bmd_fit3.mat')));
fit4 = cell2mat(struct2cell(load('bmd_fit4.mat')));
fit5 = cell2mat(struct2cell(load('bmd_fit5.mat')));
fit6 = cell2mat(struct2cell(load('bmd_fit6.mat')));
grid = cell2mat(struct2cell(load('bmd_grid.mat')));
grid = reshape(grid, [1,50]);

[x,y] = meshgrid(grid);
c = x.*y;


surf(x,y,fit1,c);
surf(x,y,fit2,c);
surf(x,y,fit3,c);
surf(x,y,fit4,c);
surf(x,y,fit5,c);
surf(x,y,fit6,c);