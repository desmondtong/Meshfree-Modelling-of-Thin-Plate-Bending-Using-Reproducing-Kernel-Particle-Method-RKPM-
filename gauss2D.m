function [G W] = gauss2D(Nod4,n)

% Gauss points & their weights
[gi wi] = NGauss(n);

ax = Nod4(1,1); bx = Nod4(2,1);
ay = Nod4(2,2); by = Nod4(3,2); %%to find distance for gauss point

% mapping the points
G1 = gi*ones(1,n);
G2 = G1';
G3 = [G1(1:end)' G2(1:end)'];
G(:,1) = (bx-ax)/2.*G3(:,1)+(ax+bx)/2;
G(:,2) = (by-ay)/2.*G3(:,2)+(ay+by)/2; %%n sets of coordinate

% scale the weights
W1 = wi*wi'; %%wi*wj(wi')
W2 = W1(1:end); %%arrange str8 down
W = (bx-ax)/2*(by-ay)/2.*W2'; %%n weight

%plot3(G(:,1),G(:,2),W*10,'*'); axis('equal')
