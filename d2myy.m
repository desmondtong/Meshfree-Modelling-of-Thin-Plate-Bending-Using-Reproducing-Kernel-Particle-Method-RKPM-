function [d2m0yy,d2m1yy,d2m11yy,d2m2yy,d2m22yy,d2m12yy] = d2myy(d2m0yy,d2m1yy,d2m11yy,d2m2yy,d2m22yy,d2m12yy,x,xA,y,yA,dc,a,W,dWdy,d2Wdy2)
d2m0yy=d2m0yy+(dc*d2Wdy2)/a;
d2m1yy=d2m1yy-(dc*(x - xA)*d2Wdy2)/a^2;
d2m11yy=d2m11yy+(dc*(x - xA)^2*d2Wdy2)/a^3;
d2m2yy=d2m2yy- (2*dc*dWdy)/a^2 - (dc*(y - yA)*d2Wdy2)/a^2;
d2m22yy=d2m22yy+(2*dc*W)/a^3 + (dc*(y - yA)^2*d2Wdy2)/a^3 + (2*dc*(2*y - 2*yA)*dWdy)/a^3;
d2m12yy=d2m12yy+(2*dc*dWdy*(x - xA))/a^3 + (dc*(x - xA)*(y - yA)*d2Wdy2)/a^3;
end

