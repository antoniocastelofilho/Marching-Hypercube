 function MarchingHyperCube(n, First, Last, Division, Func, filename) 
   file  = fopen(filename,'w');
   fprintf(file,'%3d   1\n',n);
   fprintf(file,'%3d ',Division(1:n));
   fprintf(file,'\n\n');
   
   Delta = (Last-First)./Division;
   Pert  = random('Normal',0,1,2^n,n)*0.000001;
   Base  = InitGrid(n, Division);
   for g = 0:Base(n+1)-1      
      [Grid] = GridCoords(n,g,Base); 
      fprintf('Grid = ');
      fprintf('%3d ',Grid);
      fprintf('\n');
      [VertexManifold LabelVertexManifold NumberVertex LabelVertex] = GetVertexManifold(n,Grid,Pert,First,Delta,Func);
      [VertexEdgeManifold1 NumberEdgeManifold1 NumberEdge1] = GetEdgeManifold(n,NumberVertex,LabelVertexManifold); 
      ne = 0;
      VertexEdgeManifold = [];
      for k = 1:NumberEdge1
         if NumberEdgeManifold1(k) > 2
             [Edge1 Edge2] = SolveAmbiguity(n,Grid,Pert,First,Delta,Func,LabelVertexManifold,VertexEdgeManifold1(k,:)); 
             ne = ne+1; 
             VertexEdgeManifold = [VertexEdgeManifold; VertexEdgeManifold1(k,Edge1(1)) VertexEdgeManifold1(k,Edge1(2))];
             ne = ne+1; 
             VertexEdgeManifold = [VertexEdgeManifold; VertexEdgeManifold1(k,Edge2(1)) VertexEdgeManifold1(k,Edge2(2))];
         else
             ne = ne+1; 
             VertexEdgeManifold = [VertexEdgeManifold; VertexEdgeManifold1(k,1:2)];
         end
      end
      NumberEdge = ne;
      if NumberVertex > 0
         fprintf(file,'%3d ',g);
         fprintf(file,'%3d ',Grid);
         fprintf(file,'\n');
         SkeletonHyperCube(file,n,NumberVertex,NumberEdge,VertexManifold,LabelVertexManifold,LabelVertex,VertexEdgeManifold);
         fprintf(file,'\n');
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


function [Vert FVert] = GenVert(n, Grid, Pert, First, Delta) 
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


function [Face m] = GenFace(n,k) 
   Face = []; 
   for i = 1:n-k
      Division(i) = 2*n;
   end
   [Base]  = InitGrid(n-k, Division);
   m = 0;
   for g = 0:Base(n-k+1)-1
      [F]      = GridCoords(n-k,g,Base);
      [isface] = IsFace(n,n-k,F);
      if isface == 1
         Face = [Face; F];
         m    = m+1;
      end
   end
   return;
end


function  [VertexManifold LabelVertexManifold numvert LVertex] = GetVertexManifold(n, Grid, Pert, First, Delta, Func)
   VertexManifold = [];
   LabelVertexManifold = [];
   LVertex = [];
   numvert = 0;  
   for i = 1:n-1
      Division(i) = 2*n;
   end
   [Base]  = InitGrid(n-1, Division);
   for g = 0:Base(n)-1
      [F]      = GridCoords(n-1,g,Base);
      [isface] = IsFace(n,n-1,F);
      if isface == 1
         [Label] = LabelVertex(n,F);
         [V1 FV1] = GenVertHyperCube(n,Grid,Pert,First,Delta,Label(1),Func); 
         [V2 FV2] = GenVertHyperCube(n,Grid,Pert,First,Delta,Label(2),Func); 
         A = [1 1; FV1 FV2];      
         b = [1; 0];
         lamb = A\b;
         if lamb > 0
            V = lamb(1)*V1+lamb(2)*V2;           
            VertexManifold = [VertexManifold; V];
            LabelVertexManifold = [LabelVertexManifold; F];
            LVertex = [LVertex; Label];
            numvert = numvert+1;
         end               
      end
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


function  [Label] = LabelVertex(n,edge)
   sumcoord = 0;
   for i = 1:n
      sumcoord = sumcoord + i-1;
   end
   label = 0;
   sumexp = 0;
   for i = 1:n-1
      exp    = mod(edge(i),n);
      coord  = fix(edge(i)/n);
      label  = label + coord*2^exp;
      sumexp = sumexp + exp;
   end
   Label(1) = label;
   Label(2) = label + 2^(sumcoord-sumexp);
   return
end


function [isface] = IsFace(dim,n,F)
   for i = 1:n-1
      for j = i+1:n
         if F(i) >= F(j)
            isface = 0;
            return;        
         end
         if mod(F(i),dim) == mod(F(j),dim)
            isface = 0;
            return;        
         end
      end
   end 
   isface = 1;
   return;
end

function  [VertexEdgeManifold NumberEdgeManifold NumberEdge] = GetEdgeManifold(n,NumberVertex,LabelVertexManifold)  
   NumberEdgeManifold = []; 
   VertexEdgeManifold = [];
   NumberEdge         = 0;
   for i = 1:n-1
      Division(i) = 2*n;
   end
   [Base]  = InitGrid(n-2, Division);
   for g = 0:Base(n-1)-1
      [F]      = GridCoords(n-2,g,Base);
      [isface] = IsFace(n,n-2,F);
      if isface == 1
         [LabelVertex NEdge] = InFace(n,NumberVertex,LabelVertexManifold,F);
         if NEdge == 2
             LabelVertex = [LabelVertex -1 -1];
         end
         if NEdge > 0    
            NumberEdgeManifold = [NumberEdgeManifold; NEdge];
            VertexEdgeManifold = [VertexEdgeManifold; LabelVertex];
            NumberEdge         = NumberEdge+1;
         end 
      end
   end    
   return
end

function [LabelVertex NumberEdge] = InFace(n,NumberVertex,LabelVertexManifold,F)
  LabelVertex = [];
  NumberEdge  = 0;
  for i = 1:NumberVertex     
      [k] = Include(n-2,F,n-1,LabelVertexManifold(i,:));
      if k == n-2
          LabelVertex = [LabelVertex i];
          NumberEdge = NumberEdge+1;
      end
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


function  [Edge1 Edge2] = SolveAmbiguity(n, Grid, Pert, First, Delta, Func, LabelVertexManifold, VertexEdgeManifold)
   [Label1] = LabelVertex(n,LabelVertexManifold(VertexEdgeManifold(1),:));
   [Label2] = LabelVertex(n,LabelVertexManifold(VertexEdgeManifold(2),:));
   [Label3] = LabelVertex(n,LabelVertexManifold(VertexEdgeManifold(3),:));
   [Label4] = LabelVertex(n,LabelVertexManifold(VertexEdgeManifold(4),:));                    
   Edges = [Label2; Label3; Label4];              
   [V1 FV1] = GenVertHyperCube(n,Grid,Pert,First,Delta,Label1(1),Func); 
   if (FV1 > 0) 
       [I J] = find(Edges == Label1(1));
       Edge1 = [1 I(1)+1];
       if Edge1(2) == 2
           Edge2 = [3 4];
       else
           if Edge1(2) == 3
               Edge2 = [2 4];
           else
               Edge2 = [2 3];
           end
       end
   else
       [I J] = find(Edges == Label1(2));
       Edge1 = [1 I(1)+1];
       if Edge1(2) == 2
           Edge2 = [3 4];
       else
           if Edge1(2) == 3
               Edge2 = [2 4];
           else
               Edge2 = [2 3];
           end
       end
   end      
   return
end


function SkeletonHyperCube(file, n, nv, ne, Vertex, FaceVertex, LabelVertex, Edges)
    A     = zeros(nv);
    for i = 1:ne
        A(Edges(i,1),Edges(i,2)) = 1;
        A(Edges(i,2),Edges(i,1)) = 1;
    end
    [nc nvc vc] = Split_Edges(nv,A);
    fprintf(file,'%3d\n',nc);
    for k = 1:nc
        nVertex = nvc(k);
        fprintf(file,'%3d\n',nVertex);
        for j = 1:nVertex
            FaceVertexComp(j,:) = FaceVertex(vc(k,j),:);
            LabelVertexComp(j,:) = LabelVertex(vc(k,j),:);
            VertexComp(j,:)     = Vertex(vc(k,j),:);
            fprintf(file,'%3d ',LabelVertexComp(j,:));
            fprintf(file,'%15.8f ',VertexComp(j,:));
            fprintf(file,'\n');
        end
        [nSkel Skel AdjSkel] = SkeletonComp(n, nVertex, FaceVertexComp);
        for j = 2:n
            fprintf(file,'%3d\n',nSkel(j));
            for i = 1:nSkel(j)
                fprintf(file,'%3d ',AdjSkel{j}{i}(:));
                fprintf(file,'\n');
            end
        end    
    end
    return
end


function [nc nvc vc] = Split_Edges(nv,A)
    nc   = 0;
    [k]  = FindNextComp(nv,A);
    while k ~= 0
        nc             = nc+1;
        nvc(nc)        = 1;
        vc(nc,nvc(nc)) = k;
        av             = 1;
        while av <= nvc(nc)
            k = vc(nc,av);
            for i = 1:nv
                if A(k,i) ~= 0
                    [value] = InComp(nvc(nc),i,vc(nc,:));
                    if value == 0
                        nvc(nc)        = nvc(nc)+1;
                        vc(nc,nvc(nc)) = i;
                    end
                    A(k,i) = 0;
                    A(i,k) = 0;
                end
            end
            av = av+1;
        end
        [k]  = FindNextComp(nv,A);
    end
    return
end


function [k] = FindNextComp(n, A)
    for i = 1:n-1
        for j = i+1:n
            if A(i,j) ~= 0
                k = i;
                return
            end
        end
    end
    k = 0;
    return
end


function [value] = InComp(n, y, x)
    value = 0;
    for i = 1:n
        if x(i) == y
            value = 1;
            return
        end
    end
    return
end


function [nSkel Skel AdjSkel] = SkeletonComp(n, nVertex, FaceVertex)   
    nSkel(1) = nVertex;
    Skel{1}  = FaceVertex;
    for k = 2:n-1
        for i = 1:n-k
            Division(i) = 2*n;
        end
        [Base]    = InitGrid(n-k, Division);
        FacesSimp = [];
        FacesAdj  = [];
        nSkel(k)  = 0;
        for g = 0:Base(n-k+1)-1
            [Face]    = GridCoords(n-k,g,Base);
            [isface]  = IsFace(n,n-k,Face);
            if isface == 1
                cont = 0;
                Adj  = [];
                for i = 1:nSkel(k-1)
                    [np] = Include(n-k,Face,n-k+1,Skel{k-1}(i,:));
                    if np == n-k
                        cont      = cont+1;
                        Adj(cont) = i;
                    end
                end
                if cont > 0              
                    nSkel(k)           = nSkel(k)+1;
                    FacesSimp          = [FacesSimp; Face];
                    FacesAdj{nSkel(k)} = Adj;
                end    
            end
        end
        Skel{k}    = FacesSimp;
        AdjSkel{k} = FacesAdj;
    end 
    nSkel(n)      = 1;
    AdjSkel{n}{1} = [1:nSkel(n-1)];
    return
end