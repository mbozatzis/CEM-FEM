clear
clc

a = 0.76;
b = 1.75;
V = 1;

% Geometry Description
gd = [1 1; 0 0; 0 0; b a];
ns = [82 82; 49 50];
dl = decsg(gd, 'R1-R2', ns);

% Create Mesh
[p, e, t] = initmesh(dl);
[p, e, t] = refinemesh(dl, p, e, t);
[p, e, t] = refinemesh(dl, p, e, t);
[p, e, t] = refinemesh(dl, p, e, t);
Num_nodes = size(p, 2);
Num_elements = size(t, 2);


% Define known and unknown potentials
node_id = ones(Num_nodes, 1);
Potentials = zeros(Num_nodes, 1);
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
        
        if (radius(1) < (a+b)/2)
            Potentials(nodes(1)) = V;
        else
            Potentials(nodes(1)) = 0;
        end
        if (radius(2) < (a+b)/2)
            Potentials(nodes(2)) = V;
        else
            Potentials(nodes(2)) = 0;
        end
    end
end

% Indexes of unknowns
non_pot_index = zeros(Num_nodes, 1);
counter = 0;
for in = 1:Num_nodes
    if (node_id(in) == 1)
        counter = counter + 1;
        non_pot_index(in) = counter;
    end
end

% Matrix Calculation
Nf = counter;
Sff = spalloc(Nf, Nf, 7*Nf);
B = zeros(Nf, 1);

for ie = 1:Num_elements
    nodes(1:3) = t(1:3, ie);
    region = t(4, ie);
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

            if (node_id(nodes(i)) == 1)
                if (node_id(nodes(j)) == 1)
                    Sff(non_pot_index(nodes(i)), non_pot_index(nodes(j))) = ...
                    Sff(non_pot_index(nodes(i)), non_pot_index(nodes(j))) + Se(i, j);
                else
                    B(non_pot_index(nodes(i))) = B(non_pot_index(nodes(i))) - ...
                    Se(i, j)*Potentials(nodes(j));
                end
            end
        end
    end
end

% Unknown Potential Calculation 
pot = Sff\B;
for in = 1:Num_nodes
    if node_id(in) == 1
        Potentials(in) = pot(non_pot_index(in));
    end
end

% Calculate Electric Field
[Ex, Ey] = pdegrad(p, t, Potentials);
Ex = -Ex;
Ey = -Ey;

% Plot the Results
figure(1);
pdeplot(p, e, t, 'xydata', Potentials);
axis tight;
axis equal;
hold on;
pdeplot(p, e, t, 'FlowData', [Ex', Ey']');
colormap jet;

