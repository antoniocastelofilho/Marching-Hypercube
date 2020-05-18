function MarchingSimplex(n, k, First, Last, Division, Func, filename) 
   file  = fopen(filename,'w');
   fprintf(file,'%3d %3d\n',n,k);
   fprintf(file,'%3d ',Division(1:n));
   fprintf(file,'\n\n');

   Delta = (Last-First)./Division;
   Pert  = random('Normal',0,1,2^n,n)*0.000001;
   Base  = InitGrid(n, Division);
   
   for i = 1:n
        DivisionP(i) = i;
   end
   [BaseP] = InitGrid(n,DivisionP);
   
   for i = 1:k+1
        DivisionC(i) = n-k+i;
   end
   [BaseC] = InitGrid(k+1,DivisionC);
   
   Ncube    = 0;
   Nsimplex = 0;
   for g = 0:Base(n+1)-1 
       
      Ncube  = Ncube + 1; 
      [Grid] = GridCoords(n,g,Base);
      [Vert FVert] = GenVert(n, Grid, Pert, First, Delta, Func);
        
      for s = 0:BaseP(n+1)-1
          
        Nsimplex = Nsimplex + 1;  
        fprintf('Simplexo: Ncube = %d  Nsimplex = %d  g = %d  s = %d\n',Ncube,Nsimplex,g,s);
          
        [f] = GridCoords(n,s,BaseP);
        f = f+1;
        [P] = Map_Perm(n,f);
        [Simp] = GenLabelSimplex(n,P);
        
        VertexManifold = [];
        FaceVertex     = [];
        NumberVertex   = 0;
           
        for j = 0:BaseC(k+2)-1
            
          [C]   = GridCoords(k+1,j,BaseC);
          C     = C+1;
          [lex] = Lexico(k+1,C);
          
          if lex == 1
              for i = 1:k+1
                  Face(i) = Simp(C(i));
              end 
              [Vertex trans] = GetVertexManifold(n,k,Face,Vert,FVert);
              if trans == 1
                 VertexManifold = [VertexManifold; Vertex];
                 FaceVertex     = [FaceVertex; Face];
                 NumberVertex = NumberVertex + 1;
              end
          end
          
        end
        
        if NumberVertex > 0
            
           [nSkel Skel AdjSkel] = Skeleton(n,Simp,k,NumberVertex,FaceVertex);  
           
           fprintf(file,'%3d ',g);
           fprintf(file,'%3d ',Grid);
           fprintf(file,'\n');
           fprintf(file,'  1\n');

           fprintf(file,'%3d \n',nSkel(1));
           for i = 1:nSkel(1)
              fprintf(file,'%3d ',Skel{1}(i,:));
              fprintf(file,'%15.8f ',VertexManifold(i,:));
              fprintf(file,'\n');
           end
           for j = 1:n-k
              fprintf(file,'%3d\n',nSkel(j+1));
              for i = 1:nSkel(j+1)
                 fprintf(file,'%3d ',AdjSkel{j}{i}(:));
                 fprintf(file,'\n');
              end
           end
           fprintf(file,'\n');
           
         end
        
      end
                   
   end 
   fprintf(file,'-1\n');
   
   fclose(file);
   
   return
end


function [Base] = InitGrid(n, Division) 
   Base(1) = 1;
   for i = 2:n+1
      Base(i) = Base(i-1)*Division(i-1);
   end
   return
end

    
function [Grid] = GridCoords(n,i,Base) 
   copy = i;
   for j = n:-1:2
      aux     = mod(copy,Base(j));
      Grid(j) = (copy-aux)/Base(j);
      copy    = aux;
   end
   Grid(1) = copy;
   return
end  


function [Vert FVert] = GenVert(n, Grid, Pert, First, Delta, Func) 
   Vert  = [];
   FVert = [];
   for i = 0:2^n-1 
      [Coords]    = HyperCubeCoords(n,i);
      [CoordPert] = HyperCubePert(n,Grid,Coords,Pert);
      [VHC]       = HyperCube(n,First,Delta,Grid,Coords,CoordPert);
      Vert        = [Vert; VHC];
      [FVHC]      = Func(n,VHC); 
      FVert       = [FVert; FVHC];
   end
   return
end


function [Coords] = HyperCubeCoords(n,i) 
   copy = i;
   for j = 1:n-1
      Coords(j) = mod(copy,2);
      copy      = (copy-Coords(j))/2;
   end
   Coords(n) = copy;
   return
end


function [CoordPert] = HyperCubePert(n,Grid,Coords,Pert)
   pot = 1;
   label = 0;
   for i = 1:n
      p     = abs(Coords(i) - mod(Grid(i),2));
      label = label + pot*p;
      pot   = 2*pot;    
   end
   CoordPert = Pert(label+1,:);
   return
end


function [VHC] = HyperCube(n,First,Delta,I,Coords,Pert) 
   for i = 1:n
      VHC(i) = First(i) + (I(i)+Coords(i))*Delta(i) + Pert(i);
   end
   return
end


function [Vertex trans] = GetVertexManifold(n,k,Face,Vert,FVert)        
   for i = 1:k+1
      A(1,i) = 1;
      A(2:k+1,i) = FVert(Face(i)+1,:);
   end
   b = [1; zeros(k,1)];
   lamb = A\b;
   trans = 0;
   Vertex = zeros(1,n);
   if lamb >= 0      
       for i = 1:k+1
          Vertex = Vertex + lamb(i)*Vert(Face(i)+1,:);
       end
       trans = 1;
   end   
   return
end


function [V FV] = GenVertHyperCube(n,Grid,Pert,First,Delta,Label,Func)     
   [Coords]    = HyperCubeCoords(n,Label);
   [CoordPert] = HyperCubePert(n,Grid,Coords,Pert);
   [V]         = HyperCube(n,First,Delta,Grid,Coords,CoordPert);
   [FV]        = Func(n,V); 
   return
end


function [F] = Map_Perm(n,f)
    for i = 1:n
        F(i) = 0;
    end
    for i = 1:n
        if F(f(i)) == 0
            F(f(i)) = i;
        else
            for j = n:-1:f(i)+1
                F(j) = F(j-1);
            end
            F(f(i)) = i;
        end
    end
    return
end

function [lex] = Lexico(k,f)
    for i = 1:k-1
        if f(i) >= f(i+1)
            lex = 0;
            return
        end
    end
    lex = 1;
    return
end

function [Simp] = GenLabelSimplex(n,P)
    Simp(1) = 0;
    for i = 1:n
        Simp(i+1) = Simp(i)+2^(P(i)-1);
    end
end

function [nSkel Skel AdjSkel] = Skeleton(n, Simp, k, nVertex, FaceVertex)   
    nSkel(1) = nVertex;
    Skel{1}  = FaceVertex;
    for s = k+1:n
        for i = 1:s+1
            DivisionC(i) = n-s+i;
        end
        [BaseC]      = InitGrid(s+1,DivisionC);
        FacesSimp    = [];
        FacesAdj     = [];
        nSkel(s-k+1) = 0;
        for j = 0:BaseC(s+2)-1   
          [C]   = GridCoords(s+1,j,BaseC);
          C     = C+1;
          [lex] = Lexico(s+1,C);
          if lex == 1
              for i = 1:s+1
                  Face(i) = Simp(C(i));
              end 
              cont = 0;
              Adj  = [];
              for i = 1:nSkel(s-k)
                  [np] = Include(s,Skel{s-k}(i,:),s+1,Face);
                  if np == s
                      cont      = cont+1;
                      Adj(cont) = i;
                  end
              end
              if cont > 0              
                 nSkel(s-k+1)           = nSkel(s-k+1)+1;
                 FacesSimp              = [FacesSimp; Face];
                 FacesAdj{nSkel(s-k+1)} = Adj;
              end    
           end
        end
        Skel{s-k+1}  = FacesSimp;
        AdjSkel{s-k} = FacesAdj;
    end
    return
end
              

function [k] = Include(n,x,m,y)
   k = 0;
   for i = 1:n
       j = i;
       while (j < m) && (y(j) < x(i))
           j = j+1;
       end
       if y(j) == x(i)
           k = k+1;
       end
   end
   return
end




