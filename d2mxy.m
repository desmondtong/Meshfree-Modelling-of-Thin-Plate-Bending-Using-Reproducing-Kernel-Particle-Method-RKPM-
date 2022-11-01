function [d2m0xy,d2m1xy,d2m11xy,d2m2xy,d2m22xy,d2m12xy] = d2mxy(d2m0xy,d2m1xy,d2m11xy,d2m2xy,d2m22xy,d2m12xy,x,xA,y,yA,dc,a,W,dWdx,dWdy,d2Wdxdy)
d2m0xy=d2m0xy+(dc*d2Wdxdy)/a;
d2m1xy=d2m1xy- (dc*dWdy)/a^2 - (dc*(x - xA)*d2Wdxdy)/a^2;
d2m11xy=d2m11xy+(dc*(x - xA)^2*d2Wdxdy)/a^3 + (dc*(2*x - 2*xA)*dWdy)/a^3;
d2m2xy=d2m2xy- (dc*dWdx)/a^2 - (dc*(y - yA)*d2Wdxdy)/a^2;
d2m22xy=d2m22xy+(dc*(y - yA)^2*d2Wdxdy)/a^3 + (dc*(2*y - 2*yA)*dWdx)/a^3;
d2m12xy=d2m12xy+(dc*W)/a^3 + (dc*dWdx*(x - xA))/a^3 + (dc*dWdy*(y - yA))/a^3 + (dc*(x - xA)*(y - yA)*d2Wdxdy)/a^3;
end

