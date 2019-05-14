function [ f ] = AssembleVector( ti )
%ASSEMBLEVECTOR Assembles the element vector using the current time value
global p t s np nt delta
f = zeros(np,1);

for i = 1:nt
    %Create element vectors
    fe = abs(delta(i))/6 * [s(p(t(i,1), :), ti); s(p(t(i,2), :), ti); s(p(t(i,3), :), ti)];
    %Assembly of f
    f(t(i,:)) = f(t(i,:)) + fe; 
end

end

