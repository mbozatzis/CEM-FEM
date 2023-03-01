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

% Define known and unknown fields
node_id = ones(Num_nodes, 1);
Ez = zeros(Num_nodes, 1);
Num_edges = size(e, 2);
for ie = 1:Num_edges
    region1 = e(6, ie);
    region2 = e(7, ie);
    nodes(1:2) = e(1:2, ie);
    x(1:2) = p(1, nodes(1:2));
    y(1:2) = p(2, nodes(1:2));
    radius(1:2) = sqrt(x(1:2).^2 + y(1:2).^2);

    if (region1 == 0 || region2 == 0) 
        node_id(nodes(1)) = 0;
        node_id(nodes(2)) = 0;
        Ez(nodes(1)) = 0;
        Ez(nodes(2)) = 0;
    end
end

% Indexes of unknowns
index = zeros(Num_nodes, 1);
counter = 0;
for in = 1:Num_nodes
    if (node_id(in) == 1)
        counter = counter + 1;
        index(in) = counter;
    end
end

% Matrix Calculation
Nf = counter;
S = spalloc(Nf, Nf, 7*Nf);
T = spalloc(Nf, Nf, 7*Nf);

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

            if(node_id(nodes(i)) == 1)
                if(node_id(nodes(j)) == 1)
                    S(index(nodes(i)), index(nodes(j))) = S(index(nodes(i)), index(nodes(j))) + Se(i, j);
                    T(index(nodes(i)), index(nodes(j))) = T(index(nodes(i)), index(nodes(j))) + Te(i, j);
                end
            end
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

    for in = 1:Num_nodes
        if node_id(in) == 1
            Ez(in) = V(index(in), i);
        end
    end

    figure(i-1);
    pdeplot(p, e, t, 'xydata', Ez);
    title('fc = ', fc(i-1));
    axis tight;
    axis equal;
    hold on;
    colormap jet;
end