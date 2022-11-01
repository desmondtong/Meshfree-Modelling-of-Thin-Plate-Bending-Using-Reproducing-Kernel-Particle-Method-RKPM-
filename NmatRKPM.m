function [Nmat,dNmatdx,dNmatdy,d2Nmatdx2,d2Nmatdy2,d2Nmatdxdy] = NmatRKPM(n,n1,XYi,XYui,a,dc)
Nmat=zeros(n,n1);%column must same with row size of ui,row is number of coordinate to interpolate
for i=1:n
    x=XYi(i,1);%point of interest
    y=XYi(i,2);
    m0=0;m1=0;m11=0;m2=0;m22=0;m12=0;
    dm0x=0;dm1x=0;dm11x=0;dm2x=0;dm22x=0;dm12x=0;dm0y=0;dm1y=0;dm11y=0;dm2y=0;dm22y=0;dm12y=0;
    d2m0xx=0;d2m1xx=0;d2m11xx=0;d2m2xx=0;d2m22xx=0;d2m12xx=0;
    d2m0yy=0;d2m1yy=0;d2m11yy=0;d2m2yy=0;d2m22yy=0;d2m12yy=0;
    d2m0xy=0;d2m1xy=0;d2m11xy=0;d2m2xy=0;d2m22xy=0;d2m12xy=0;
    for j=1:n1
        xA=XYui(j,1);%existing data
        yA=XYui(j,2);
        [W,dWdx,dWdy,d2Wdx2,d2Wdy2,d2Wdxdy] = wfuncRKPM_2D(xA,x,yA,y,a);
        m0=m0+(1/a)*W*dc;
        m1=m1+((xA-x)/a)*(1/a)*W*dc;
        m11=m11+((xA-x)/a)^2*(1/a)*W*dc;
        m2=m2+((yA-y)/a)*(1/a)*W*dc;
        m22=m22+((yA-y)/a)^2*(1/a)*W*dc;
        m12=m12+((yA-y)/a)*((xA-x)/a)*(1/a)*W*dc;
        
        %1st derivative
        [dm0x,dm1x,dm11x,dm2x,dm22x,dm12x] = dmx(dm0x,dm1x,dm11x,dm2x,dm22x,dm12x,x,xA,y,yA,dc,a,W,dWdx);
        [dm0y,dm1y,dm11y,dm2y,dm22y,dm12y] = dmy(dm0y,dm1y,dm11y,dm2y,dm22y,dm12y,x,xA,y,yA,dc,a,W,dWdy);
        
        %2nd derivative
        [d2m0xx,d2m1xx,d2m11xx,d2m2xx,d2m22xx,d2m12xx] = d2mxx(d2m0xx,d2m1xx,d2m11xx,d2m2xx,d2m22xx,d2m12xx,x,xA,y,yA,dc,a,W,dWdx,d2Wdx2);
        [d2m0yy,d2m1yy,d2m11yy,d2m2yy,d2m22yy,d2m12yy] = d2myy(d2m0yy,d2m1yy,d2m11yy,d2m2yy,d2m22yy,d2m12yy,x,xA,y,yA,dc,a,W,dWdy,d2Wdy2);
        [d2m0xy,d2m1xy,d2m11xy,d2m2xy,d2m22xy,d2m12xy] = d2mxy(d2m0xy,d2m1xy,d2m11xy,d2m2xy,d2m22xy,d2m12xy,x,xA,y,yA,dc,a,W,dWdx,dWdy,d2Wdxdy);
    end
    
    D=m0*(m11*m22-m12^2)-(m1^2*m22-2*m1*m2*m12+m2^2*m11);
    C1=(m11*m22-m12^2)/D;
    C21=(m1*m22-m2*m12)/D;
    C22=(m2*m11-m1*m12)/D;
    
    %1st derivative
    [dC1x,dC21x,dC22x] = dCx(m0,m1,m11,m2,m22,m12,dm0x,dm1x,dm11x,dm2x,dm22x,dm12x);
    [dC1y,dC21y,dC22y] = dCy(m0,m1,m11,m2,m22,m12,dm0y,dm1y,dm11y,dm2y,dm22y,dm12y);
    
    %2nd derivative
    [d2C1xx,d2C21xx,d2C22xx] = d2Cxx(m0,m1,m11,m2,m22,m12,dm0x,dm1x,dm11x,dm2x,dm22x,dm12x,d2m0xx,d2m1xx,d2m11xx,d2m2xx,d2m22xx,d2m12xx);
    [d2C1yy,d2C21yy,d2C22yy] = d2Cyy(m0,m1,m11,m2,m22,m12,dm0y,dm1y,dm11y,dm2y,dm22y,dm12y,d2m0yy,d2m1yy,d2m11yy,d2m2yy,d2m22yy,d2m12yy);
    [d2C1xy,d2C21xy,d2C22xy] = d2Cxy(m0,m1,m11,m2,m22,m12,dm0x,dm1x,dm11x,dm2x,dm22x,dm12x,dm0y,dm1y,dm11y,dm2y,dm22y,dm12y,d2m0xy,d2m1xy,d2m11xy,d2m2xy,d2m22xy,d2m12xy);
    
    for k=1:n1
        xA=XYui(k,1);
        yA=XYui(k,2);
        [W,dWdx,dWdy,d2Wdx2,d2Wdy2,d2Wdxdy] = wfuncRKPM_2D(xA,x,yA,y,a);
        Exs=C1+C21*((x-xA)/a)+C22*((y-yA)/a);
        Nmat(i,k)=Exs*(1/a)*W*dc;
        
        %1st derivative
        dExsdx=dC1x + C21/a + (dC21x*(x - xA))/a + (dC22x*(y - yA))/a;
        dExsdy=dC1y + C22/a + (dC21y*(x - xA))/a + (dC22y*(y - yA))/a;
        
        dNmatdx(i,k)=dExsdx*(1/a)*W*dc+Exs*(1/a)*dWdx*dc; %NR eq55a Liu1995
        dNmatdy(i,k)=dExsdy*(1/a)*W*dc+Exs*(1/a)*dWdy*dc;
        
        %2nd derivative
        d2Exsdx2=(2*dC21x)/a + ((x - xA)*d2C21xx)/a + ((y - yA)*d2C22xx)/a + d2C1xx;
        d2Exsdy2=(2*dC22y)/a + ((x - xA)*d2C21yy)/a + ((y - yA)*d2C22yy)/a + d2C1yy;
        d2Exsdxdy=dC22x/a + dC21y/a + ((x - xA)*d2C21xy)/a + ((y - yA)*d2C22xy)/a + d2C1xy;
        
        d2Nmatdx2(i,k)=(dc*Exs*d2Wdx2)/a + (dc*W*d2Exsdx2)/a + (2*dc*dExsdx*dWdx)/a;
        d2Nmatdy2(i,k)=(dc*Exs*d2Wdy2)/a + (dc*W*d2Exsdy2)/a + (2*dc*dExsdy*dWdy)/a;
        d2Nmatdxdy(i,k)=(dc*Exs*d2Wdxdy)/a + (dc*W*d2Exsdxdy)/a + (dc*dExsdx*dWdy)/a + (dc*dExsdy*dWdx)/a;
    end
end
end

