function curlyE = expectation_functions(X,Pi,a_i,a_pi,T)
%----------------------------------------------------------------
% Description: Given a function x(e,a) of the state, calclutates
%        expectation functions {eps1,...,epsT} up to horizon T.
%              eps{j}(e,a)=E[ eps{j-1}(e',a')|(e,a) ]
%           where eps{0}(e,a)=x(e,a).
%-----------------------------------------------------------------
% Set up array of curlyEs and fill in first row with base case
    curlyE=zeros([T,size(X)]);
    curlyE(1,:,:) = X;
% Recursively apply law of iterated expectations
    for j = 2:T
        tempX=squeeze(curlyE(j-1,:,:));
        curlyE(j,:,:) = expectation_iteration(tempX,Pi,a_i,a_pi);
    end
end

function X_out = expectation_iteration(X, Pi, a_i, a_pi)
% Matrix multiplication using *
 Xend = Pi * X;
% Call the expectation_policy function to get the result
 X_out = expectation_policy(Xend, a_i, a_pi);
end

function X = expectation_policy(Xend, a_i, a_pi)
    % 先获取a_i的维度
    [m, n] = size(a_i);
    % 初始化X为与Xend相同大小的零矩阵
    X = zeros(size(Xend));
    % 双重循环遍历a_i的每一个元素
    for e = 1:m
        for a = 1:n
            % 考虑MATLAB索引从1开始
            index=a_i(e,a)+1;
            % 如果index或index+1超出Xend的列数，则需要特殊处理
            if index < size(Xend,2)
                X(e,a)= a_pi(e,a)*Xend(e,index) + (1-a_pi(e,a))*Xend(e,index+1);
            elseif index >= size(Xend,2)
                X(e,a)= a_pi(e,a)*Xend(e,index);
            end
        end
    end
end