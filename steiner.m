function Inew = steiner(I, a, m)
    Inew = I;
    for i = 1:3
        for j= 1:3
            Inew(i,j) = Inew(i,j) - m*a(i)*a(j);
            if i==j
                Inew(i,j) = Inew(i,j) + m*norm(a)^2;
            end
        end
    end
end
