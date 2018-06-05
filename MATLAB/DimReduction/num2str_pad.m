function y = num2str_pad(x)


x_str = num2str(x);

k = length( x_str );

if k == 1
    y = ['00' x_str ];
elseif k == 2
    y = ['0' x_str];
elseif k == 3
    y = x_str;
end

end


