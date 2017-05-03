function zad2(I, N)
    syms a b l3 t real
    m = sym('m%d', [1, 3]);
    q = sym('q%d', [3, 1]);

    x1 = [0; 0; 0];
    x2 = [0; 0; b];
    x3 = [-a*sin(q(1))+q(3)*cos(q(2))*cos(q(1)); a*cos(q(1))+q(3)*cos(q(2))*sin(q(1)); b-q(3)*sin(q(2))];
    fi1 = [0; 0; q(1)];
    fi2 = [-q(2)*sin(q(1)); q(2)*cos(q(1)); q(1)];
    fi3 = fi2;

    J = sym(zeros(3, 3, 3));
    J(3, 3, 1) = sym('I1z');
    J(1, 1, 2) = sym('I2');
    J(2, 2, 2) = sym('I2');
    J(3, 3, 2) = sym('I2');
    J(2, 2, 3) = sym('I3y');
    J(3, 3, 3) = sym('I3z');

    k = 0.4 * I;
    n = 0.5 * N; % zmienilem z "m" na "n", zeby nie kolidowalo z masa
    c = 0.15 * I;
    d = 0.3 * N;
    e = 1.5 * I;
    f = 0.6 * N;
    h = 0.5 * N;

    bval = 0.15;
    m1 = 1.5;
    J1z = 0.215;

    aval = 0.09;
    m2 = 1.8;
    J2 = 0.35;

    m3 = 2;
    l3val = 0.45;

    qval = [k*pi/6 + cos(n*t^2); c*pi/4 - cos(d*t^3); 0.01*e - 0.03*cos(f*t^2) - l3val];
    tval = h;
    mval = [m1 m2 m3];

    Jval = zeros(3, 3, 3);
    Jval(3, 3, 1) = J1z;
    Jval(1, 1, 2) = J2;
    Jval(2, 2, 2) = J2;
    Jval(3, 3, 2) = J2;
    Jval(2, 2, 3) = (m3 * l3val^2) / 12;
    Jval(3, 3, 3) = (m3 * l3val^2) / 12;
    % nie licze Jval(1, 1, 3), bo podobno belka pryzmatyczna ma znikomy ten trzeci moment bezwladnosci

    common(m, J, q, x1, x2, x3, fi1, fi2, fi3, qval, Jval, [t m l3 a], [tval mval l3val aval]);
end
