function zad2(I, N)
    syms a b l3 t J1z J2xyz J3yz J3z positive
    m = sym('m%d', [1, 3], 'positive');
    q = sym('q%d', [3, 1], 'real');

    x1 = [0; 0; 0];
    x2 = [0; 0; b];
    x3 = [-a*sin(q(1))+q(3)*cos(q(2))*cos(q(1)); a*cos(q(1))+q(3)*cos(q(2))*sin(q(1)); b-q(3)*sin(q(2))];
    fi1 = [0; 0; q(1)];
    fi2 = [-q(2)*sin(q(1)); q(2)*cos(q(1)); q(1)];
    fi3 = fi2;

    J = sym(zeros(3, 3, 3));
    J(3, 3, 1) = J1z;
    
    J2 = sym(zeros(3, 3));
    J2(1, 1) = J2xyz;
    J2(2, 2) = J2xyz;
    J2(3, 3) = J2xyz;
    R01 = [cos(q(1)), -sin(q(1)), 0; sin(q(1)), cos(q(1)), 0; 0, 0, 1];
    J(:, :, 2) = R01 * J2 * R01';
    
    J3(2, 2) = J3yz;
    J3(3, 3) = J3yz;
    % pomijam J3(1, 1), bo podobno belka pryzmatyczna ma znikomy ten trzeci moment bezwladnosci
    R02 = R01 * [cos(q(2)), 0, sin(q(2)); 0, 1, 0; -sin(q(2)), 0, cos(q(2))];
    J(:, :, 3) = R02 * J3 * R02';

    k = 0.4 * I;
    n = 0.5 * N; % zmienilem z "m" na "n", zeby nie kolidowalo z masa
    c = 0.15 * I;
    d = 0.3 * N;
    e = 1.5 * I;
    f = 0.6 * N;
    h = 0.5 * N;

    m1 = 1.5;
    J1zval = 0.215;

    aval = 0.09;
    m2 = 1.8;
    J2xyzval = 0.35;

    m3 = 2;
    l3val = 0.4; % wybralem losowo ten z dolu strony

    qval = [k*pi/6 + cos(n*t^2); c*pi/4 - cos(d*t^3); 0.01*e - 0.03*cos(f*t^2) - l3val];
    tval = h;
    mval = [m1 m2 m3];

    J3yzval = (m3 * l3val^2) / 12;

    common(m, J, q, x1, x2, x3, fi1, fi2, fi3, qval, [t m l3 a J1z J2xyz J3yz], [tval mval l3val aval J1zval J2xyzval J3yzval]);
end
