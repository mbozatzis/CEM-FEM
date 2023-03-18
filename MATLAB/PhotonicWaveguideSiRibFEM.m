clear
clc

% TE mode Calculation

% distanses in um (10^(-6))
w = 0.4 * 10^(-6); 
H = 0.34 * 10^(-6); 
h = 0.05 * 10^(-6); 
A = 4;% * 10^(-6);
R = 1;% * 10^(-6);
e1 = 12;
e2 = 3.9;
m0 = 1.256637061*10^(-6);
e0 = 8.854187817* 10^(-12); 
c0 = 3*10^8;


% Geometry Description 
rect1 = [3; 4; -A/2; -A/2; A/2; A/2; -A/2; A/2; A/2; -A/2];
C2 = [1; 0; 0; R; 0; 0; 0; 0; 0; 0];
rect2 = [3; 4; -A/2; -A/2; A/2; A/2; -A/16; -A/16 + h; -A/16 + h; -A/16];
rect3 = [3; 4; -w/2; -w/2; w/2; w/2; -A/16 + h; -A/16 + h + H; -A/16 + h + H; -A/16 + h];
rect4 = [3; 4; -A/2; -A/2; A/2; A/2; -A/2; -A/16; -A/16; -A/2];

gd = [rect1, rect2, rect3, rect4];
gd = [rect1, C2];
sf = 'rect1 + rect2 + rect3 + rect4';
sf = 'rect1 + C2';
ns = char('rect1', 'rect2', 'rect3', 'rect4');
ns = char('rect1', 'C2');
ns = ns';

% Create Mesh
dl = decsg(gd, sf, ns);
[p,e,t] = initmesh(dl); 
[p,e,t] = refinemesh(dl,p,e,t); 
[p,e,t] = refinemesh(dl,p,e,t);
Num_nodes = size(p, 2);
Num_elements = size(t, 2);


% Matrix Calculation
Sxx = spalloc(Num_nodes, Num_nodes, 7*Num_nodes);
Sxy = spalloc(Num_nodes, Num_nodes, 7*Num_nodes);
Sxz = spalloc(Num_nodes, Num_nodes, 7*Num_nodes);
Syx = spalloc(Num_nodes, Num_nodes, 7*Num_nodes);
Syy = spalloc(Num_nodes, Num_nodes, 7*Num_nodes);
Syz = spalloc(Num_nodes, Num_nodes, 7*Num_nodes);
Szx = spalloc(Num_nodes, Num_nodes, 7*Num_nodes);
Szy = spalloc(Num_nodes, Num_nodes, 7*Num_nodes);
Szz = spalloc(Num_nodes, Num_nodes, 7*Num_nodes);
Txx = spalloc(Num_nodes, Num_nodes, 7*Num_nodes);
Tyy = spalloc(Num_nodes, Num_nodes, 7*Num_nodes);
Tzz = spalloc(Num_nodes, Num_nodes, 7*Num_nodes);

for ie = 1:Num_elements
    nodes(1:3) = t(1:3, ie);
    x(1:3) = p(1, nodes(1:3));
    y(1:3) = p(2, nodes(1:3));
    area = t(4, ie);

    De = det([1 x(1) y(1); 1 x(2) y(2); 1 x(3) y(3)]);
    Ae = abs(De/2);
    a(1) = (x(2)*y(3) - x(3)*y(2))/De;
    a(2) = (x(3)*y(1) - x(1)*y(3))/De;
    a(3) = (x(1)*y(2) - x(2)*y(1))/De;
    b(1) = (y(2) - y(3))/De;
    b(2) = (y(3) - y(1))/De;
    b(3) = (y(1) - y(2))/De;
    c(1) = (x(3) - x(2))/De;
    c(2) = (x(1) - x(3))/De;
    c(3) = (x(2) - x(1))/De;

    for i = 1:3
        for j = 1:3
            if (area == 2 || area == 3)
                er = e1;
            elseif area == 4
                er = e2;
            else
                er = 1;
            end


            Sxye(i, j) = -(1/sqrt(er*e0/m0))*c(i)*b(j)*Ae;
            Sxze(i, j) = -1i*(c(j)*1/3)*Ae;
            Syxe(i, j) = -(1/sqrt(er*e0/m0))*b(i)*c(j)*Ae;
            Syze(i, j) = -1i*(b(j)*1/3)*Ae;
            Szxe(i, j) = -1i*(c(i)*1/3)*Ae;
            Szye(i, j) = -1i*(b(i)*1/3)*Ae;
            Szze(i, j) = (1/sqrt(er*e0/m0))*(2*c(i)*c(j))*Ae;

            if i == j
                Te(i, j) = (1/sqrt(er*e0*m0))*er*e0*Ae/6 + sqrt(er*e0/m0)*Ae/6;
                Tz(i, j) = (1/sqrt(er*e0*m0))*er*e0*Ae/6;
                Sxxe(i, j) = (-sqrt(er*e0/m0)*(1/(6))+ (1/sqrt(er*e0/m0))*c(i)*c(j))*Ae;
                Syye(i, j) = (-sqrt(er*e0/m0)*(1/(6))+ (1/sqrt(er*e0/m0))*b(i)*b(j))*Ae;
            else
                Te(i, j) = (1/sqrt(er*e0*m0))*er*e0*Ae/12 + sqrt(er*e0/m0)*Ae/12;
                Tz(i, j) = (1/sqrt(er*e0*m0))*er*e0*Ae/12;
                Sxxe(i, j) = (-sqrt(er*e0/m0)*(1/(12)+ (1/sqrt(er*e0/m0))*c(i)*c(j)))*Ae;
                Syye(i, j) = (-sqrt(er*e0/m0)*(1/(12)+ (1/sqrt(er*e0/m0))*b(i)*b(j)))*Ae;
            end
             Sxxe(i, j) = (1/sqrt(er*e0/m0))*(c(i)*c(j))*Ae;
             Syye(i, j) = (1/sqrt(er*e0/m0))*(b(i)*b(j))*Ae;

            Sxx(nodes(i), nodes(j)) = Sxx(nodes(i), nodes(j)) + Sxxe(i, j);
            Sxy(nodes(i), nodes(j)) = Sxy(nodes(i), nodes(j)) + Sxye(i, j);
            Sxz(nodes(i), nodes(j)) = Sxz(nodes(i), nodes(j)) + Sxze(i, j);
            Syx(nodes(i), nodes(j)) = Syx(nodes(i), nodes(j)) + Syxe(i, j);
            Syy(nodes(i), nodes(j)) = Syy(nodes(i), nodes(j)) + Syye(i, j);
            Syz(nodes(i), nodes(j)) = Syz(nodes(i), nodes(j)) + Syze(i, j);
            Szx(nodes(i), nodes(j)) = Szx(nodes(i), nodes(j)) + Szxe(i, j);
            Szy(nodes(i), nodes(j)) = Szy(nodes(i), nodes(j)) + Szye(i, j);
            Szz(nodes(i), nodes(j)) = Szz(nodes(i), nodes(j)) + Szze(i, j);
            Txx(nodes(i), nodes(j)) = Txx(nodes(i), nodes(j)) + Te(i, j);
            Tyy(nodes(i), nodes(j)) = Tyy(nodes(i), nodes(j)) + Te(i, j);
            Tzz(nodes(i), nodes(j)) = Tzz(nodes(i), nodes(j)) + Tz(i, j);
        end
    end
end

S = [Sxx, Sxy, Sxz; Syx, Syy, Syz; Szx, Szy, Szz];
T = [Txx, zeros(size(Txx)), zeros(size(Txx)); zeros(size(Txx)), Tyy, zeros(size(Txx)); ...
    zeros(size(Txx)), zeros(size(Txx)), Tzz];

% Field and cut-off frequency calculation
k = 6;
sigma = 'smallestabs';
[V, D] = eigs(S, T, k, sigma);

pos = length(V)/3;
Ex = real(V(1:pos, :));
Ey = real(V(pos+1:2*pos, :));
Ez = imag(V(2*pos+1:3*pos, :));

for i = 1:6
    subplot(2, 3, i);
    pdeplot(p, e, t, 'xydata', Ex(:, i));
    axis tight;
    axis equal;
    hold on;
    colormap jet;
end

