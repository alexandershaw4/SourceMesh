function y=vec2meshvals(vi,vals)
%
%
%

y = vi*0;
for i = 1:length(vals)
    f = find(vi==i);
    y(f) = vals(i);
end