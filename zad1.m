function zad1(I, N)
    syms lc3 lc2 l2 t J2z J3x J3y J3z real
    m = sym('m%d', [1, 3], 'real');
    q = sym('q%d', [3, 1], 'real');

    x1 = [q(1); 0; 0];
    x2 = [q(1); 0; -lc2];
    x3 = [q(1)+lc3*sin(q(3))*cos(q(2)); lc3*sin(q(3))*sin(q(2)); -l2-lc3*cos(q(3))];
    fi1 = [0; 0; 0];
    fi2 = [0; 0; q(2)];
    fi3 = [q(3)*sin(q(2)); -q(3)*cos(q(2)); q(2)];

    J = sym(zeros(3, 3, 3));
    J(3, 3, 2) = J2z;
    R = [cos(q(2)), -sin(q(2)), 0; sin(q(2)), cos(q(2)), 0; 0, 0, 1];
    J3 = sym(zeros(3,3));
    J(1, 1) = J3x;
    J(2, 2) = J3y;
    J(3, 3) = J3z;
    J(:, :, 3) = R * J3 * R';

    a = 0.01 * I;
    b = 0.02 * N;
    c = 0.5 * I;
    d = 0.025 * N;
    e = 0.05 * I;
    f = 0.06 * N;
    h = 0.5 * N;

    m1 = 1;

    m2 = 1.2;
    J2zval = 0.025;

    m31 = 0.4;
    r31 = 0.16 / 2;
    l31 = 0.02;
    m32 = 1.8;
    l32 = 0.6;
    lc3val = l32 * m32 / (2 * (m31 + m32)); % wzor na srodek ciezkosci
    display(lc3val)

    qval = [a + b*t^2; c*pi/6 + d*t^3; e*pi/12 + f*t^2];
    tval = h;
    mval = [m1 m2 m31+m32];

    Jval31 = zeros(3);
    Jval31(1, 1) = (3*r31^2 + l31^2) * m31 / 12;
    Jval31(2, 2) = (m31 * r31^2) / 2;
    Jval31(3, 3) = (3*r31^2 + l31^2) * m31 / 12;
    Jval32 = zeros(3);
    Jval32(1, 1) = (m32 * l32^2) / 12;
    Jval32(2, 2) = (m32 * l32^2) / 12;
    % nie licze Jval32(3,3), bo podobno belka pryzmatyczna ma znikomy ten trzeci moment bezwladnosci

    J3val = steiner(Jval31, [0, 0, -lc3val], m31) + steiner(Jval32, [0, 0, l32/2 - lc3val], m32);
    % przesuwam obydwa tensory do srodka ciezkosci calego czlonu trzeciego

    common(m, J, q, x1, x2, x3, fi1, fi2, fi3, qval, [t m lc3 J2z reshape(J3, [1, 9])], [tval mval lc3val J2zval reshape(J3val, [1, 9])]);
end
