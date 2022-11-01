function [dm0y,dm1y,dm11y,dm2y,dm22y,dm12y] = dmy(dm0y,dm1y,dm11y,dm2y,dm22y,dm12y,x,xA,y,yA,dc,a,W,dWdy)
dm0y=dm0y+dc*dWdy/a;
dm1y=dm1y-(dc*dWdy*(x - xA))/a^2;
dm11y=dm11y+(dc*dWdy*(x - xA)^2)/a^3;
dm2y=dm2y-(dc*W)/a^2 - (dc*dWdy*(y - yA))/a^2;
dm22y=dm22y+(dc*dWdy*(y - yA)^2)/a^3 + (dc*W*(2*y - 2*yA))/a^3;
dm12y=dm12y+(dc*W*(x - xA))/a^3 + (dc*dWdy*(x - xA)*(y - yA))/a^3;
end

