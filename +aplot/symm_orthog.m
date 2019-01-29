function y = symm_orthog(x)
% Efficient symmetric orthogonalisation method for matrices, using:
%
% y = x * real(inv(x' * x)^(1/2));
%
% AS

fprintf('\nOrthogonalising\n');
y = [x] * real(inv([x]' * [x])^(1/2));
y = (max(x(:)) - min(x(:))) * ( (y - min(y(:))) / (max(y(:)) - min(y(:))) );

end