function y = row(x)
% make vector row if it is not
if iscolumn(x)
    y = x.';
else
    y = x;
end
