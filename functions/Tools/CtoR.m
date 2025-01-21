function Y = CtoR(X)
% Convert complex-valued input to real output while keeping matrix orientation

if ~isreal(X) % Only convert if input is complex
    % Alsways put largest dimension as rows
    if size(X,2) > size(X,1)
        X = X.';
        transposeFlag = 1;
    else
        transposeFlag = 0;
    end

    % Convert to complex numbers
    Y = zeros(size(X,1),size(X,2)*2);
    for n = 1:size(X,2)
        Y(:,1+(n-1)*2:2+(n-1)*2) = [real(X(:,n)) imag(X(:,n))];
    end

    % Convert back to original orientation
    if transposeFlag
        Y = Y.';
    end
else
    Y = X;
end

end