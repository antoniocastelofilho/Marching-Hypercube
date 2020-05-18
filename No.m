function No
  
    First    = [-1.7 -1.7 -1.1];
    Last     = [1.7 1.7 1.1];
    Division = [28 28 22]; 
  
    %MarchingSimplex(3, 1, First, Last, Division, @no_p, 'no_30_MS.pol');
    MarchingHyperCube(3, First, Last, Division, @no_p, 'no_30_MHC.pol');

    
    return 

    
    function [f] = ESFERA(x,x0,r) 
        f  = (x(1)-x0(1))*(x(1)-x0(1))+(x(2)-x0(2))*(x(2)-x0(2))+(x(3)-x0(3))*(x(3)-x0(3))-r*r;
        return
    end

    function [x0] = PARAMETRIZACAO(t) 
       x0(1) = (1.0+0.4*cos(1.5*t))*cos(t);
       x0(2) = (1.0+0.4*cos(1.5*t))*sin(t);
       x0(3) = 0.5*sin(1.5*t);
   
%      x0(1) = (8.0+3.0*cos(5.0*t))*cos(2.0*t);
%      x0(2) = (8.0+3.0*cos(5.0*t))*sin(2.0*t);
%      x0(3) = 3.0*sin(5.0*t);
       return
    end

    function [f] = no_p(n,x) 
      p    = 30.0;
      k  = 150;
      r  = 0.24;
      f(1) = 0.0;
      for i = 0:k
         t    = 4.0*pi*i/k;
         [x0] = PARAMETRIZACAO(t);    
         f(1) = f(1) + (ESFERA(x,x0,r)+1.0)^(-p);
      end
      f(1) = (f(1))^(-1.0/p)-1.0;
    end
    
 end 
