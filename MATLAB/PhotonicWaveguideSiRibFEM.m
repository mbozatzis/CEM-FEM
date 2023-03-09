clear
clc

% TE mode Calculation

% distanses in um (10^(-6))
w = 0.4; 
H = 0.34; 
h = 0.05; 
A = 4;
B = 4;
e1 = 12;
e2 = 3.9;

% Geometry Description 
rect1 = [3; 4; -A/2; -A/2; A/2; A/2; -A/2; A/2; A/2; -A/2];
rect2 = [3; 4; -A/2; -A/2; A/2; A/2; -A/16; -A/16 + h; -A/16 + h; -A/16];
rect3 = [3; 4; -w/2; -w/2; w/2; w/2; -A/16 + h; -A/16 + h + H; -A/16 + h + H; -A/16 + h];
rect4 = [3; 4; -A/2; -A/2; A/2; A/2; -A/2; -A/16; -A/16; -A/2];

gd = [rect1, rect2, rect3, rect4];
sf = 'rect1 + rect2 + rect3 + rect4';
ns = char('rect1', 'rect2', 'rect3', 'rect4');
ns = ns';

% Create Mesh
dl = decsg(gd, sf, ns);
[p,e,t] = initmesh(dl); 
[p,e,t] = refinemesh(dl,p,e,t); 
Num_nodes = size(p, 2);
Num_elements = size(t, 2);

pdeplot(p,e,t);



