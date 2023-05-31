
function indices = compute_reduced_indices(before, after)
% Compute the retained functional (colour) values after patch reduction
%

indices = zeros(length(after), 1);
for i = 1:length(after)
    dotprods = (before * after(i, :)') ./ sqrt(sum(before.^2, 2));
    [~, indices(i)] = max(dotprods);
end
end