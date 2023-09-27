function polyMatrix = c_gen_PolyMatrix(data, maxDegree)
% Description:
% 输入 data       包含多列不同变量的矩阵
%      maxDegree  多项式最大次数
% 输出 polyMatrix 包含交乘项的多项式矩阵，每一列代表一项

% Example:
% A,B,C 3列,2次的情况下:
%    A B C A^2 AB AC B^2 BC C^2
%------------------------------------------------------
   [numRows, numCols] = size(data);
    
    % 初始化输出矩阵为输入数据矩阵
    polyMatrix = data;

    % 定义一个递归函数来生成笛卡尔积
    function addColumns(currentDegree, currentProduct, colIndex)
        if currentDegree == 0
            if isempty(currentProduct)
                return;
            end
            [~, loc] = ismember(currentProduct', polyMatrix', 'rows');
            if loc == 0
                polyMatrix = [polyMatrix, currentProduct];
            end
            return;
        end
        
        for i = colIndex:numCols
            newProduct = currentProduct .* data(:, i);
            addColumns(currentDegree - 1, newProduct, i);
        end
    end

    for degree = 2:maxDegree
        addColumns(degree, ones(numRows, 1), 1);
    end
end