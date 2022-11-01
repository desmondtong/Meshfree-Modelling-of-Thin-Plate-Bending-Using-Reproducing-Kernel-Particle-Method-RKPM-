clear, clc, clf

% define the structures properties
L = 120; %mm
B = 1.0*L; %ratio b/a
tw = 1; %mm
BPts = [0 0; L 0; L B; 0 B]; % boundary points

v = 0.3; % poissons ratio
E = 210000; %N/mm2
D = (tw^3/12)*E/(1-v^2)*[1 v 0; v 1 0; 0 0 (1-v)/2];

% Support conditions [restraints]
alpha = 1e8*E;
SuPts = [0 0 , L 0 , alpha 0 , 0 -1;
    L 0 , L B , alpha 0 , 1 0;
    L B , 0 B , alpha 0 , 0 1;
    0 B , 0 0 , alpha 0 , -1 0]; % [node1 node2 a1 a2 n1 n2]

% Applied Loading
pz = -1; % pressure (kN/m^2)

% Variable parameters
div = 4;
h1 = L/div;
h2 = B/round(B/h1); %square mesh
a = h1;

clear xi
% Generate nodes in the cantilever
[xi1, xi2] = meshgrid(0:h1:L,0:h2:B);
xi(:,1) = xi1(1:end);
xi(:,2) = xi2(1:end); %%exisiting point
% plot(xi(:,1),xi(:,2),'gp'); hold on; axis('equal');


%% Background cells
%to get gauss weight and abiscisa
Nc1 = round(1.0*L/h1);
Nc2 = round(1.0*B/h2); %%1.0 to control size of grid (div)
clear xc
[xi1, xi2] = meshgrid(0:L/Nc1:L,0:B/Nc2:B);
xc = [xi1(:),xi2(:)]; %%point of interest

Bcells = zeros(4,2,1); k = 0; %%do layer..but actually can combine str8 down to btm like xQ %%one layer one coordinate of grid
for i = 1:Nc1
    for j = 1:Nc2
        k = k + 1;
        l = j+(Nc2+1)*(i-1);
        C = [l l+Nc2+1 l+Nc2+2 l+1];
        Bcells(:,:,k) = xc(C,:);
        %TT = xc(C,:);
        %plot(TT([1:end 1],1),TT([1:end 1],2),'-*'); hold on
    end
end

Ng = 2; %%no of Gauss point
xQ=[];
WxQ=[]; %%open kosong one/when cat() str8 replace
for l = 1:Nc1*Nc2
    % loop for cells here
    [G, W] = gauss2D(Bcells(:,:,l),Ng); %Gauss points
    xQ=cat(1,xQ,G);
    WxQ=cat(1,WxQ,W); %%grouping all G and W in one vector
    %plot(G(:,1),G(:,2),'r.'); hold on
end
%plot(xQ(:,1),xQ(:,2),'ko');


%% RKPM shape function for all xQ
[Nmat,dNmatdx,dNmatdy,d2Nmatdx2,d2Nmatdy2,d2Nmatdxdy] = NmatRKPM(length(xQ),length(xc),xQ,xc,a,h1);


%% Stiffness Matrix Kij

clear f K
f = zeros(size(xi,1),1);
K = zeros(size(xi,1));
for nQ = 1:size(xQ,1)
    f = f + WxQ(nQ)*tw*Nmat(nQ,:)'*pz;
    Bj = [-d2Nmatdx2(nQ,:)' -d2Nmatdy2(nQ,:)' -2*d2Nmatdxdy(nQ,:)']';
    K = K + WxQ(nQ)*Bj'*D*Bj; %%actually need double loop for bi bj, use vector faster
end


%% Penalty Matrix Ka
Ka = zeros(size(xi,1));
n = 9; %no of gauss point
h = 1; %division of boundaries
for Nsu = 1:size(SuPts,1) %%integrate over boundary only not domain
    % loop over edges/boundary (4)
    xQs=zeros(n*h,2);
    if Nsu == 1 || Nsu == 3
        b1=SuPts(Nsu,1);
        b2=SuPts(Nsu,3);
        [xQs1,WxQs] = gauss_1D(h,n,b1,b2);
        xQs(:,1)=xQs1;
        xQs(:,2)=SuPts(Nsu,2);
    else
        b1=SuPts(Nsu,2);
        b2=SuPts(Nsu,4);
        [xQs1,WxQs] = gauss_1D(h,n,b1,b2);
        xQs(:,2)=xQs1;
        xQs(:,1)=SuPts(Nsu,1);
    end
    Alpha = diag([SuPts(Nsu,5);SuPts(Nsu,6)]);
    [Nmats,dNmatdxs,dNmatdys] = NmatRKPM(length(xQs),length(xc),xQs,xc,a,h1);
    for nQ = 1:size(xQs,1)
        % loop over the quadrature (gauss) points
        Ns = [Nmats(nQ,:)' zeros(length(xi),1)]';
        Ka = Ka + WxQs(nQ)*Ns'*Alpha*Ns;
    end
end
K = K + Ka;


%% wertical displacements wnj
U = K\f;
fprintf('f = %.5f\n',sum(f)) %pz*A=f
fprintf('A = %.5f\n',sum(WxQ))
fprintf('u = %.5f\n',max(abs(U)))
fprintf('a = %.5f\n',max(abs(U*D(1)/(pz*L^4)))) %alpha in timoshenko
fprintf('p = %.5f\n',max(abs(D(1)/(pz*L^4))))


%% plot results
Scale = 1;

plot3(0,0,0) %initial plot
hold on
plot3(xi(:,1),xi(:,2),zeros(size(xi(:,2))),'b+')
plot3(xi(:,1),xi(:,2),U/Scale,'r.')
plot(BPts([1:end 1],1),BPts([1:end 1],2),'b-')
axis equal
hold off

