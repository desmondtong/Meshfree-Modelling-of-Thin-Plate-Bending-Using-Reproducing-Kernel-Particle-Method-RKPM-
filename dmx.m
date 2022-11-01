function [dm0x,dm1x,dm11x,dm2x,dm22x,dm12x] = dmx(dm0x,dm1x,dm11x,dm2x,dm22x,dm12x,x,xA,y,yA,dc,a,W,dWdx)
dm0x=dm0x+dWdx*dc/a;
dm1x=dm1x-(dc*W)/a^2-(dc*dWdx*(x-xA))/a^2;
dm11x=dm11x+(dc*dWdx*(x-xA)^2)/a^3+(dc*W*(2*x-2*xA))/a^3;
dm2x=dm2x-(dc*dWdx*(y - yA))/a^2;
dm22x=dm22x+(dc*dWdx*(y - yA)^2)/a^3;
dm12x=dm12x+(dc*W*(y - yA))/a^3 + (dc*dWdx*(x - xA)*(y - yA))/a^3;
end

