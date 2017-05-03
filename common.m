function common(m, J, q, x1, x2, x3, fi1, fi2, fi3, qval, oldval, newval)
    syms g qd1 qd2 qd3 qdd1 qdd2 qdd3 real
    qd = sym('qd%d', [3, 1], 'real');
    qdd = sym('qdd%d', [3, 1], 'real');
    
    Jv1 = jacobian(x1, q);
    display(Jv1)
    Jv2 = jacobian(x2, q);
    display(Jv2)
    Jv3 = jacobian(x3, q);
    display(Jv3)
    Jo1 = jacobian(fi1, q);
    display(Jo1)
    Jo2 = jacobian(fi2, q);
    display(Jo2)
    Jo3 = jacobian(fi3, q);
    display(Jo3)
    
    D = m(1) * (Jv1'*Jv1) + m(2) * (Jv2'*Jv2) + m(3) * (Jv3'*Jv3) +  Jo1'*J(:,:,1)*Jo1 + Jo2'*J(:,:,2)*Jo2 + Jo3'*J(:,:,3)*Jo3;

    C = sym(zeros(3, 3, 3));
    for i = 1:3
        for j = 1:3
            for k = 1:3
                C(i, j, k) = (diff(D(j, k), q(i)) + diff(D(i, k), q(j)) + diff(D(i, j), q(k))) / 2;
            end
        end
    end
    
    V = g * (m(1) * x1(3) + m(2) * x2(3) + m(3) * x3(3));
    
    f = gradient(V, q);
    for k = 1:3
        for j = 1:3
            f(k) = f(k) + D(k, j) * qdd(j);
            for i = 1:3
                f(k) = f(k) + C(i, j, k) * qd(i) * qd(j);
            end
        end
    end
    
    D = simplify(D);
    C = simplify(C);
    V = simplify(V);
    f = simplify(f, 'Steps', 200);
    display(D)
    display(C)
    display(V)
    display(f)

    qdval = diff(qval);
    qddval = diff(qdval);
    gval = 9.81;
    display(double(subs(subs(f, [q qd qdd], [qval qdval qddval]), [oldval, g], [newval, gval])))
end
