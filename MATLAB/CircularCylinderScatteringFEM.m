clear
clc

lambda = 1;
f = 300*10^6;
k = (2*pi)/lambda;
a = lambda/8;
R = a + lambda/2;
E0 = 1;
m0 = 1.256637061*10^(-6);
e0 = 8.854187817* 10^(-12);  

alpha = -(1/(2*R) + 1i*k);
omega = ((2*pi*f)^2);

% Geometry Description
gd = [1 1; 0 0; 0 0; R a];
ns = [82 82; 49 50];
dl = decsg(gd, 'R1-R2', ns);

% Create Mesh
[p, e, t] = initmesh(dl);
[p, e, t] = refinemesh(dl, p, e, t);
[p, e, t] = refinemesh(dl, p, e, t);
[p, e, t] = refinemesh(dl, p, e, t);
Num_nodes = size(p, 2);
Num_elements = size(t, 2);

pdeplot(p, e, t);

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
        
        if (radius(1) < R/2) 
            Ez(nodes(1)) = -E0*exp(-1i*2*pi*p(1,nodes(1)));
            node_id(nodes(1)) = 0;
        else
            node_id(nodes(1)) = 2;
        end
        if (radius(2) < R/2) 
            Ez(nodes(2)) = -E0*exp(-1i*2*pi*p(1,nodes(2)));
            node_id(nodes(2)) = 0;
        else
            node_id(nodes(2)) = 2;
        end
    end
end

% Indexes of unknowns
index = zeros(Num_nodes, 1);
counter = 0;
for in = 1:Num_nodes
    if (node_id(in) ~= 0)
        counter = counter + 1;
        index(in) = counter;
    end
end


% Matrix Calculation
Nf = counter;
Sff = spalloc(Nf, Nf, 7*Nf);
Tff = spalloc(Nf, Nf, 7*Nf);
Tffc = spalloc(Nf, Nf, ceil(2*pi*sqrt(Nf)));
B = zeros(Nf, 1);


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
            Se(i, j) = (1/m0)*(b(i)*b(j) + c(i)*c(j))*Ae;
            if i == j
                Te(i, j) = e0*Ae/6;
            else
                Te(i, j) = e0*Ae/12;
            end

            if(node_id(nodes(i)) == 1)
                if(node_id(nodes(j)) ~= 0)
                    Sff(index(nodes(i)), index(nodes(j))) = Sff(index(nodes(i)), index(nodes(j))) + Se(i, j);
                    Tff(index(nodes(i)), index(nodes(j))) = Tff(index(nodes(i)), index(nodes(j))) - omega*Te(i, j);
                else
                    B(index(nodes(i))) = B(index(nodes(i))) - (Se(i, j) - omega*Te(i, j))*Ez(nodes(j));
                end
            end
            if(node_id(nodes(i)) == 2)
                if(node_id(nodes(j)) == 1)
                    Sff(index(nodes(i)), index(nodes(j))) = Sff(index(nodes(i)), index(nodes(j))) + Se(i, j);
                    Tff(index(nodes(i)), index(nodes(j))) = Tff(index(nodes(i)), index(nodes(j))) - omega*Te(i, j);
                elseif(node_id(nodes(j)) == 2)
                    Sff(index(nodes(i)), index(nodes(j))) = Sff(index(nodes(i)), index(nodes(j))) + Se(i, j);
                    Tff(index(nodes(i)), index(nodes(j))) = Tff(index(nodes(i)), index(nodes(j))) - omega*Te(i, j);
                    Tffc(index(nodes(i)), index(nodes(j))) = Tffc(index(nodes(i)), index(nodes(j))) + (1/(m0*e0))*alpha*Te(i,j);
                end
            end
        end
    end
end

A = Sff + Tff + Tffc;
el = A\B;
for in = 1:Num_nodes
    if (node_id(in) ~= 0)
        Ez(in) = abs(el(index(in))+exp(-1i*k*p(1,in))); 
     else
         Ez(in) = abs(Ez(in)+exp(-1i*k*p(1,in)));
    end
end

pdeplot(p, e, t, 'xydata', Ez);
axis tight;
axis equal;
hold on;
colormap jet;
