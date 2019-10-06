%clc;
clear all;
close all;

function e = element(p1, p2)
    s = p2 - p1;
    n = [s(2), -s(1)];
    e = struct('p1', p1, 'p2', p2, 'pc', 0.5 * (p1 + p2), 'n', n / norm(n), 'u', 0, 'q', 0);
end

function g = green(p1, p2)
    r = p2 - p1;
    g = -log(norm(r)) / (2 * pi);
end

function gg = greenGradient(p1, p2)
    r = p2 - p1;
    gg = - r / (dot(r, r) * 2 * pi);
end

function showElements(els)
    N = length(els);
    f1 = figure(1);
    for i = 1 : N
        plot([els{i}.p1(1) els{i}.p2(1)], [els{i}.p1(2) els{i}.p2(2)], 'k', 'linewidth', 2);
        if i == 1
            hold on;
        end
        plot(els{i}.pc(1), els{i}.pc(2), 'ro');
        plot([els{i}.pc(1) els{i}.pc(1) + 0.1 * els{i}.n(1)], [els{i}.pc(2) els{i}.pc(2) + 0.1 * els{i}.n(2)], 'b', 'linewidth', 2);
    end
    axis([-2, 2 -2, 2], 'square');
end

N = 80;
A = 2;
B = 1;
phi = linspace(2 * pi, 0, N + 1);
x = A * cos(phi);
y = B * sin(phi);

elems_phi1 = cell(N, 1);
for i = 1 : N
    elems_phi1{i} = element([x(i) y(i)], [x(i + 1) y(i + 1)]);
end
%showElements(elems_phi1);

elems_phi2 = elems_phi1;
elems_phi3 = elems_phi1;

for i = 1 : N
    elems_phi1{i}.u = -elems_phi1{i}.n(1);
    elems_phi2{i}.u = -elems_phi2{i}.n(2);
    elems_phi3{i}.u = -(elems_phi3{i}.pc(1) * elems_phi3{i}.n(2) - elems_phi3{i}.pc(2) * elems_phi3{i}.n(1));
end

H = zeros(N);
G = zeros(N);
phi1 = zeros(N, 1);
dphi1 = zeros(N, 1);
phi2 = zeros(N, 1);
dphi2 = zeros(N, 1);
phi3 = zeros(N, 1);
dphi3 = zeros(N, 1);
max_l = 0;
for i = 1 : N
    phi1(i) = elems_phi1{i}.u;
    phi2(i) = elems_phi2{i}.u;
    phi3(i) = elems_phi3{i}.u;
    for j = 1 : N
        if i == j
            H(i, j) = 0.5;
            r1 = 0.5 * norm(elems_phi1{i}.p2 - elems_phi1{i}.p1);
            G(i, j) = r1 * (1 - log(r1)) / pi;
        else
            s = elems_phi1{j}.p2 - elems_phi1{j}.p1;
            l = norm(s);
            s = s / l;
            dG = @(z)green(elems_phi1{i}.pc, elems_phi1{j}.p1 + z * s);
            dH = @(z)dot(greenGradient(elems_phi1{i}.pc, elems_phi1{j}.p1 + z * s), elems_phi1{j}.n);
            G(i, j) = quad(dG, 0, l);
            H(i, j) = quad(dH, 0, l);
            %l = norm(elems_phi1{i}.p2 - elems_phi1{i}.p1);
            %if max_l < l
            %    max_l = l;
            %end
            %H(i, j) = dot(greenGradient(elems_phi1{i}.pc, elems_phi1{j}.pc), elems_phi1{j}.n) * l;
            %G(i, j) = green(elems_phi1{i}.pc, elems_phi1{j}.pc) * l;
            %H(i, j) = (    dot(greenGradient(elems_phi1{i}.pc, elems_phi1{j}.p1), elems_phi1{j}.n) + ...
            %           4 * dot(greenGradient(elems_phi1{i}.pc, elems_phi1{j}.pc), elems_phi1{j}.n) + ...
            %               dot(greenGradient(elems_phi1{i}.pc, elems_phi1{j}.p2), elems_phi1{j}.n)) * l / 6;
            %G(i, j) = (    green(elems_phi1{i}.pc, elems_phi1{j}.p1) + ...
            %           4 * green(elems_phi1{i}.pc, elems_phi1{j}.pc) + ... 
            %               green(elems_phi1{i}.pc, elems_phi1{j}.p2)) * l / 6;
        end
    end
end
max_l
Z = G \ H;
dphi1 = Z * phi1;
dphi2 = Z * phi2;
dphi3 = Z * phi3;
%dphi1 = G \ (H * phi1);
%dphi2 = G \ (H * phi2);
%dphi3 = G \ (H * phi3);

lam1 = 0;
lam2 = 0;
lam3 = 0;
max_l = 0;
for i = 1 : N
    l = norm(elems_phi1{i}.p2 - elems_phi1{i}.p1);
    if max_l < l
        max_l = l;
    end
    elems_phi1{i}.q = dphi1(i);
    elems_phi2{i}.q = dphi2(i);
    elems_phi3{i}.q = dphi3(i);


    lam1 = lam1 + phi1(i) * dphi1(i) * l;
    lam2 = lam2 + phi2(i) * dphi2(i) * l;
    lam3 = lam3 + phi3(i) * dphi3(i) * l;

end

fprintf(1, '%.15g | %.15g\n', lam1, pi * B^2);
fprintf(1, '%.15g | %.15g\n', lam2, pi * A^2);
fprintf(1, '%.15g | %.15g\n', lam3, pi * (A^2 - B^2)^2 / 8.0);
max_l
