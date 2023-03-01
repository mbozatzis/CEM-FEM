clear
clc

w = 4;
h = 0.2;
d = 1;
er = 2.2;
V = 100;

A = 2*w;
B = 2*w;

% Geometry Description
gd = [3 3 3; 
      4 4 4; 
      -A/2 w/2 w/2; 
      A/2 w/2 w/2; 
      A/2 -w/2 -w/2; 
      -A/2 -w/2 -w/2;
      -B/2 d/2 -d/2; 
      -B/2 d/2+h -d/2-h; 
      B/2 d/2+h -d/2-h; 
      B/2 d/2 -d/2];
sf = 'R1-R2-R3';
ns = [82 82 82; 49 50 51];
dl = decsg(gd, sf, ns);

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

    if (region1 == 0 || region2 == 0)
        if  (y(1) >= d/2 - h) && (y(1) <= d/2 + h) && (x(1) >= -w/2) && (x(1) <= w/2)
            Potentials(nodes(1)) = V/2;
            node_id(nodes(1)) = 0;
        elseif (y(1) >= -d/2 - h) && (y(1) <= -d/2 + h) && (x(1) >= -w/2) && (x(1) <= w/2)
            Potentials(nodes(1)) = -V/2;
            node_id(nodes(1)) = 0;
        end
        if (y(2) >= d/2 - h) && (y(2) <= d/2 + h) && (x(2) >= -w/2) && (x(2) <= w/2)
            Potentials(nodes(1)) = V/2;
            node_id(nodes(1)) = 0;
        elseif (y(2) >= -d/2 - h) && (y(2) <= -d/2 + h) && (x(2) >= -w/2) && (x(2) <= w/2)
            Potentials(nodes(1)) = -V/2;
            node_id(nodes(1)) = 0;
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
            flag(1) = (x(i) >= -w/2 && x(i) <= w/2) && (y(i) >= -d/2 && y(i) <= d/2);
            flag(2) = (x(j) >= -w/2 && x(j) <= w/2) && (y(j) >= -d/2 && y(j) <= d/2);
            if(flag(1) && flag(2))
                Se(i,j) = er*(b(i)*b(j) + c(i)*c(j))*Ae;
            else
                Se(i, j) = (b(i)*b(j) + c(i)*c(j))*Ae;
            end

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
pdeplot(p, e, t, 'xydata', Potentials, 'contour', 'on');
axis tight;
axis equal;
hold on;
pdeplot(p, e, t, 'FlowData', [Ex', Ey']');
colormap jet;

