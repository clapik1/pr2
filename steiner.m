function I = steiner(a, m)
    I = sym(zeros(3, 3));
    for i = 1:3
        for j= 1:3
            I(i,j) = -m*a(i)*a(j);
            if i==j
                I(i,j) = I(i,j) + m*norm(a)^2;
            end
        end
    end
end
