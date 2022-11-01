function [W,dWdx,dWdy,d2Wdx2,d2Wdy2,d2Wdxdy] = wfuncRKPM_2D(xui,xi,yui,yi,a)
z=(sqrt((xi-xui)^2+(yi-yui)^2))/a;
if -2<=z && z<=-1
    W=(1/6)*(2+z)^3;
    dWdx=((2*xi - 2*xui)*(2*a + ((xi - xui)^2 + (yi - yui)^2)^(1/2))^2)/(4*a^3*((xi - xui)^2 + (yi - yui)^2)^(1/2));
    dWdy=((2*yi - 2*yui)*(2*a + ((xi - xui)^2 + (yi - yui)^2)^(1/2))^2)/(4*a^3*((xi - xui)^2 + (yi - yui)^2)^(1/2));
    d2Wdx2=(4*a^2*((xi - xui)^2 + (yi - yui)^2) + 4*a*((xi - xui)^2 + (yi - yui)^2)^(3/2) + xi^2*((xi - xui)^2 + (yi - yui)^2) + xui^2*((xi - xui)^2 + (yi - yui)^2) + ((xi - xui)^2 + (yi - yui)^2)^2 - 4*a^2*xi^2 - 4*a^2*xui^2 + 8*a^2*xi*xui - 2*xi*xui*((xi - xui)^2 + (yi - yui)^2))/(2*a^3*((xi - xui)^2 + (yi - yui)^2)^(3/2));
    d2Wdy2=(4*a^2*((xi - xui)^2 + (yi - yui)^2) + 4*a*((xi - xui)^2 + (yi - yui)^2)^(3/2) + yi^2*((xi - xui)^2 + (yi - yui)^2) + yui^2*((xi - xui)^2 + (yi - yui)^2) + ((xi - xui)^2 + (yi - yui)^2)^2 - 4*a^2*yi^2 - 4*a^2*yui^2 + 8*a^2*yi*yui - 2*yi*yui*((xi - xui)^2 + (yi - yui)^2))/(2*a^3*((xi - xui)^2 + (yi - yui)^2)^(3/2));
    d2Wdxdy=((xi - xui)*(yi - yui)*(- 4*a^2 + xi^2 - 2*xi*xui + xui^2 + yi^2 - 2*yi*yui + yui^2))/(2*a^3*((xi - xui)^2 + (yi - yui)^2)^(3/2));
elseif -1<=z && z<=0
    W=2/3-z^2*(1+z/2);
    dWdx=-((4*a + 3*((xi - xui)^2 + (yi - yui)^2)^(1/2))*(xi - xui))/(2*a^3);
    dWdy=-((4*a + 3*((xi - xui)^2 + (yi - yui)^2)^(1/2))*(yi - yui))/(2*a^3);
    d2Wdx2=- (4*a + 3*((xi - xui)^2 + (yi - yui)^2)^(1/2))/(2*a^3) - (3*(2*xi - 2*xui)*(xi - xui))/(4*a^3*((xi - xui)^2 + (yi - yui)^2)^(1/2));
    d2Wdy2=- (4*a + 3*((xi - xui)^2 + (yi - yui)^2)^(1/2))/(2*a^3) - (3*(2*yi - 2*yui)*(yi - yui))/(4*a^3*((xi - xui)^2 + (yi - yui)^2)^(1/2));
    d2Wdxdy=-(3*(2*yi - 2*yui)*(xi - xui))/(4*a^3*((xi - xui)^2 + (yi - yui)^2)^(1/2));
elseif 0<=z && z<=1
    W=2/3-z^2*(1-z/2);
    dWdx=-((4*a - 3*((xi - xui)^2 + (yi - yui)^2)^(1/2))*(xi - xui))/(2*a^3);
    dWdy=-((4*a - 3*((xi - xui)^2 + (yi - yui)^2)^(1/2))*(yi - yui))/(2*a^3);
    d2Wdx2=(3*(2*xi - 2*xui)*(xi - xui))/(4*a^3*((xi - xui)^2 + (yi - yui)^2)^(1/2)) - (4*a - 3*((xi - xui)^2 + (yi - yui)^2)^(1/2))/(2*a^3);
    d2Wdy2=(3*(2*yi - 2*yui)*(yi - yui))/(4*a^3*((xi - xui)^2 + (yi - yui)^2)^(1/2)) - (4*a - 3*((xi - xui)^2 + (yi - yui)^2)^(1/2))/(2*a^3);
    d2Wdxdy=(3*(2*yi - 2*yui)*(xi - xui))/(4*a^3*((xi - xui)^2 + (yi - yui)^2)^(1/2));
elseif 1<=z && z<=2
    W=(1/6)*(2-z)^3;
    dWdx=-((2*xi - 2*xui)*(2*a - ((xi - xui)^2 + (yi - yui)^2)^(1/2))^2)/(4*a^3*((xi - xui)^2 + (yi - yui)^2)^(1/2));
    dWdy=-((2*yi - 2*yui)*(2*a - ((xi - xui)^2 + (yi - yui)^2)^(1/2))^2)/(4*a^3*((xi - xui)^2 + (yi - yui)^2)^(1/2));
    d2Wdx2=-(4*a^2*((xi - xui)^2 + (yi - yui)^2) - 4*a*((xi - xui)^2 + (yi - yui)^2)^(3/2) + xi^2*((xi - xui)^2 + (yi - yui)^2) + xui^2*((xi - xui)^2 + (yi - yui)^2) + ((xi - xui)^2 + (yi - yui)^2)^2 - 4*a^2*xi^2 - 4*a^2*xui^2 + 8*a^2*xi*xui - 2*xi*xui*((xi - xui)^2 + (yi - yui)^2))/(2*a^3*((xi - xui)^2 + (yi - yui)^2)^(3/2));
    d2Wdy2=-(4*a^2*((xi - xui)^2 + (yi - yui)^2) - 4*a*((xi - xui)^2 + (yi - yui)^2)^(3/2) + yi^2*((xi - xui)^2 + (yi - yui)^2) + yui^2*((xi - xui)^2 + (yi - yui)^2) + ((xi - xui)^2 + (yi - yui)^2)^2 - 4*a^2*yi^2 - 4*a^2*yui^2 + 8*a^2*yi*yui - 2*yi*yui*((xi - xui)^2 + (yi - yui)^2))/(2*a^3*((xi - xui)^2 + (yi - yui)^2)^(3/2));
    d2Wdxdy=-((xi - xui)*(yi - yui)*(- 4*a^2 + xi^2 - 2*xi*xui + xui^2 + yi^2 - 2*yi*yui + yui^2))/(2*a^3*((xi - xui)^2 + (yi - yui)^2)^(3/2));
else
    W=0;
    dWdx=0;
    dWdy=0;
    d2Wdx2=0;
    d2Wdy2=0;
    d2Wdxdy=0;
end
d2Wdx2(isnan(d2Wdx2))=-2;
d2Wdy2(isnan(d2Wdy2))=-2;
d2Wdxdy(isnan(d2Wdxdy))=0;
end

%cubic spline reproducing kernel formulation
%eqn 36-38 (RKPM struc. dynamic)