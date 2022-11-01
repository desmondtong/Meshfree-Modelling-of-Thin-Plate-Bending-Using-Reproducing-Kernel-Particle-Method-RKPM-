function [d2m0xx,d2m1xx,d2m11xx,d2m2xx,d2m22xx,d2m12xx] = d2mxx(d2m0xx,d2m1xx,d2m11xx,d2m2xx,d2m22xx,d2m12xx,x,xA,y,yA,dc,a,W,dWdx,d2Wdx2)
d2m0xx=d2m0xx+(dc*d2Wdx2)/a;
d2m1xx=d2m1xx-(2*dc*dWdx)/a^2 - (dc*(x - xA)*d2Wdx2)/a^2;
d2m11xx=d2m11xx+(2*dc*W)/a^3 + (dc*(x - xA)^2*d2Wdx2)/a^3 + (2*dc*(2*x - 2*xA)*dWdx)/a^3;
d2m2xx=d2m2xx-(dc*(y - yA)*d2Wdx2)/a^2;
d2m22xx=d2m22xx+(dc*(y - yA)^2*d2Wdx2)/a^3;
d2m12xx=d2m12xx+(2*dc*dWdx*(y - yA))/a^3 + (dc*(x - xA)*(y - yA)*d2Wdx2)/a^3;
end

