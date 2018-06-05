function [m,tie] = LDA_score(M,s,S,Pi,q)

tie = 0;

for i = 1:6
    delta(i) =-prod(S(i,:)) - sum( ( ( q-M(i,:) ).^2 )./s ) + 2*log(Pi(i));
end

m = find( delta == max(delta) );

if length(m) >= 2
    tie=tie+1;
    m = m(1,1);
end
