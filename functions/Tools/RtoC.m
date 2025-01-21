function Y = RtoC(X)
% Convert real-valued input to complex output while keeping matrix orientation

if isreal(X) % Only convert if input is real
    % Alsways put largest dimension as rows
    if size(X,2) > size(X,1)
        X = X.';
        transposeFlag = 1;
    else
        transposeFlag = 0;
    end

    % Check for multiple of 2 dimensions
    if mod(size(X,2),2) ~= 0
        error('Cannot convert non-multiple of 2 dimensions to complex numbers');
    end

    % Convert to complex numbers
    Y = complex(zeros(size(X,1),size(X,2)/2));
    for n = 1:size(X,2)/2
        Y(:,n) = X(:,1+(n-1)*2) + 1i*X(:,2+(n-1)*2);
    end

    % Convert back to original orientation
    if transposeFlag
        Y = Y.';
    end
else
    Y = X;
end
end