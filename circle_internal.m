%clc;
clear all;
close all;

function showBoundary(x, y)
    N = length(x);
    f1 = figure(1);
    for i = 1 : N
        j = i + 1;
        if j > N
            j -= N;
        end
        plot([x(i) x(j)], [y(i) y(j)], 'k', 'linewidth', 2);
        if i == 1
            hold on;
        end
    end
    plot(x, y, 'ro', 'markerfacecolor', 'r');
    axis('square');
end

function showElements(elems)
    N = length(elems);
    f1 = figure(1);
    for i = 1 : N
        p1 = elems{i}.p1;
        p2 = elems{i}.p2;
        plot([p1(1) p2(1)], [p1(2), p2(2)], 'k', 'linewidth', 2);
        if i == 1
            hold on;
        end
        plot(elems{i}.pc(1), elems{i}.pc(2), 'ro', 'markerfacecolor', 'r');
        plot(p1(1), p1(2), 'bo', 'markerfacecolor', 'b');
        plot(p2(1), p2(2), 'bo', 'markerfacecolor', 'b');
    end

    axis('square');
end

function u = preciseSolution(x, y)
    %u = x .* sin(2 * pi * y) + y .* cos(3.2 * pi * x);
    u = sin(x) .* cosh(y);
end

function g = preciseGradient(x, y)
    g = zeros(1, 2);
    %g(1) = sin(2 * pi * y) - 3.2 * pi * y * sin(3.2 * pi * x);
    %g(2) = 2 * pi * x * cos(2 * pi * y) + cos(3.2 * pi * x);
    g(1) = cos(x) .* cosh(y);
    g(2) = sin(x) .* sinh(y);
end

function g = green(p1, p2)
    g = log(1 / norm(p1 - p2)) / (2 * pi);
end

function gg = greenGradient(p1, p2)
    r = p2 - p1;
    gg = -r / (norm(r)^2 * 2 * pi);
end

function showSolution()
    r = linspace(0, 1, 50);
    phi = linspace(0, 2 * pi, 101);
    [rr, pp] = meshgrid(r, phi);
    x = rr .* cos(pp);
    y = rr .* sin(pp);
    u = preciseSolution(x, y);
    f2 = figure(2);
    colormap(rainbow(64));
    [c, h] = contourf(x, y, u, 64);
    set(h, 'linecolor', 'none');
    axis('square');
end

function e = elem(p1, p2, type)
    %type = 0 -> u
    %type = 1 -> q
    s = p2 - p1;
    s /= norm(s);
    e = struct('p1', p1, 'p2', p2, 'pc', 0.5 * (p1 + p2), 'n', [s(2), -s(1)], ...
               'u', 0, 'q', 0, 'type', type);
end

N = 40;
dt = 2 * pi / N;
t = linspace(0, 2 * pi - dt, N);
A = 0.2;
B = 0.1;
x = A * cos(t);
y = B * sin(t);


elems = cell(N, 1);
for i = 1 : N
    j = i + 1;
    if j > N
        j = 1;
    end
    elems{i} = elem([x(i) y(i)], [x(j) y(j)], 0);
end

for i = 1 : N
    if elems{i}.type == 0
        elems{i}.u = preciseSolution(elems{i}.pc(1), elems{i}.pc(2));
    elseif elems{i}.type == 1
        elems{i}.q = dot(preciseGradient(elems{i}.pc(1), elems{i}.pc(2)), elems{i}.n);
    else
        fprintf(stderr, 'Unknown type of element\n');
        return;
    end
end

H = zeros(N);
G = zeros(N);
for i = 1 : N
    for j = 1 : N
        if i == j
            H(i, j) = 0.5;
            s = elems{i}.p2 - elems{i}.p1;
            r1 = norm(s) * 0.5;
            G(i, j) = r1 * (1 - log(r1)) / pi;
        else
            s = elems{j}.p2 - elems{j}.p1;
            l = norm(s);
            s = s / l;
            dG = @(z)green(elems{i}.pc, elems{j}.p1 + z * s);
            dH = @(z)dot(greenGradient(elems{i}.pc, elems{j}.p1 + z * s), elems{j}.n);
            G(i, j) = quad(dG, 0, l);
            H(i, j) = quad(dH, 0, l);
            %s = elems{j}.p2 - elems{j}.p1;
            %ns = norm(s);
            %r = elems{j}.pc - elems{i}.pc;
            %nr = norm(r);
            %G(i, j) = green(elems{i}.pc, elems{j}.pc) * ns;
            %H(i, j) = dot(greenGradient(elems{i}.pc, elems{j}.pc), elems{j}.n) * ns;
        end
    end
end


U = zeros(N, 1);
Q = zeros(N, 1);
for i = 1 : N
    U(i) = elems{i}.u;
    Q(i) = dot(preciseGradient(elems{i}.pc(1), elems{i}.pc(2)), elems{i}.n);
end
%H * U - G * Q
Q1 = G \ (H * U);
max(abs(Q1 - Q))
%for i = 1 : N
%    fprintf(1, '%10g\t%10g\n', Q1(i), Q(i));
%end


%showElements(elems);
%showBoundary(x, y);
%showSolution();
