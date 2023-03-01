clear
clc

a = 1;

% Geometry Description 
gd = [1; 0; 0; a];
dl = decsg(gd);

% Create Mesh
[p, e, t] = initmesh(dl);
[p, e, t] = refinemesh(dl, p, e, t);
[p, e, t] = refinemesh(dl, p, e, t);
[p, e, t] = refinemesh(dl, p, e, t);
Num_nodes = size(p, 2);
Num_elements = size(t, 2);

% Matrix Calculation
S = spalloc(Num_nodes, Num_nodes, 7*Num_nodes);
T = spalloc(Num_nodes, Num_nodes, 7*Num_nodes);

for ie = 1:Num_elements
    nodes(1:3) = t(1:3, ie);
    x(1:3) = p(1, nodes(1:3));
    y(1:3) = p(2, nodes(1:3));

    De = det([1 x(1) y(1); 1 x(2) y(2); 1 x(3) y(3)]);
    Ae = abs(De/2);
    b(1) = (y(2) - y(3))/De;
    b(2) = (y(3) - y(1))/De;
    b(3) = (y(1) - y(2))/De;
    c(1) = (x(3) - x(2))/De;
    c(2) = (x(1) - x(3))/De;
    c(3) = (x(2) - x(1))/De;

    for i = 1:3
        for j = 1:3
            Se(i, j) = (b(i)*b(j) + c(i)*c(j))*Ae;
            if i == j
                Te(i, j) = Ae/6;
            else
                Te(i, j) = Ae/12;
            end
            S(nodes(i), nodes(j)) = S(nodes(i), nodes(j)) + Se(i, j);
            T(nodes(i), nodes(j)) = T(nodes(i), nodes(j)) + Te(i, j);
        end
    end
end

% Field and cut-off frequency calculation
k = 6;
sigma = 'smallestabs';
[V, D] = eigs(S, T, k, sigma);

c0 = 3e8;
fc = zeros(5, 1);
for i = 2:6
    fc(i-1) = c0*sqrt(D(i,i))/(2*pi);

    figure(i-1);
    pdeplot(p, e, t, 'xydata', V(:,i));
    title('fc = ', fc(i-1));
    axis tight;
    axis equal;
    hold on;
    colormap jet;
end



