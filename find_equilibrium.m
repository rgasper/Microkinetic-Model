function cutoff = find_equilibrium(A, epsilon)
%function finds the point in a set of column diff. vectors where all of
%them are below a maximum gradient value.
sz = size(A);
check = abs(A) < epsilon;
for i = 1:sz(1)
    tf = 1;
    for j = 1:sz(2)
        tf = tf*check(i,j);
    end
    if tf
        cutoff = i;
        break
    end
end
end
