function y = identify(x);


y=zeros(1,length(x));
for i =1:length(x)
    id = x(i);
    if id <= 320
                n = 1;
            elseif id>=321 && id<=434
                n = 2;
            elseif id>=435 && id<=460
                n = 3;
            elseif id>=461 && id<=505
                n = 4;
            elseif id>=506 && id<=607
                n = 5;
            elseif id>=608 && id<=729
                n = 6;
    end
    
    y(i)=n;
end