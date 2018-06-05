function y = rand_pick(a,b,n)

i=0;
y = zeros(n,1);
k=0;
while k==0
    x = randi([a,b],1,1);
    if ~any(x == y)  
        i=i+1;
        y(i) = x;
    end
        
    if i == n
        k=1;
    end
        
end

end
