function y = col(x)
% make vector column if it is not
if isrow(x)
    y = x.';
else
    y = x;
end
