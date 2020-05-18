function Toro
    First    = [-1.5 -1.5 -1.5 -1.5];
    Last     = [1.5 1.5 1.5 1.5];
    Division = [10 10 10 10]; 
    FirstPoint = [1 0 0.5 0];
  
    MarchingSimplex(4, 1, First, Last, Division, @toro_R4, 'toro_R4_New_MS.pol');
    %MarchingHyperCube(4, First, Last, Division, @toro_R4, 'toro_R4_New_MHC.pol');
    %ContinuationSimplex(4, 1, First, Last, Division, FirstPoint, @toro_R4, 'toro_R4_New_CS.pol');

    
    return 

%     function [f] = toro_R4(n,x) 
%       f(1) = (1.0-sqrt(x(1)^2+x(2)^2+x(4)^2))^2+x(3)^2-0.25;
%       return
%     end  

    function [f] = toro_R4(n,x) 
      f(1) = (1.0-sqrt(x(1)^2+x(2)^2))^2+(1.0-sqrt(x(3)^2+x(4)^2))^2-0.25;

      return
    end    
 end 