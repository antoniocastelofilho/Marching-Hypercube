function ContinuationHypercube(n, k, First, Last, Division, FirstPoint, Func, filename)
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
    
   NHcubeTrans = 0;
   HcubeTrans  = [];
   NHcubeProc  = 0;
   HcubeProc   = [];
   
   [g value] = GetFirstHypercube(n,k,First,Last,Division,FirstPoint,Func);
   
   if value < 0
       fprintf('Nenhum simplexo transversal\n');
       return
   end
   
   [NHcubeTrans HcubeTrans] = InsertList(NHcubeTrans,HcubeTrans,g);
   
   Ncube = 0;
   while NHcubeTrans > 0
        [g NHcubeTrans HcubeTrans] = GetAndRemoveList(NHcubeTrans,g,HcubeTrans);
        [NHcubeProc HcubeProc]     = InsertList(NHcubeProc,HcubeProc,g);
        Ncube = Ncube + 1;
        fprintf('Hypercube: Ncube = %d  g = %d  N = %d\n',Ncube,g,NHcubeTrans);
        [Grid]                     = GridCoords(n,g,Base);
        [Vert FVert]               = GenVert(n, Grid, Pert, First, Delta, Func);
              
        [VertexManifold LabelVertexManifold NumberVertex LabelVertex] = GetVertexManifold(n,Grid,Pert,First,Delta,Func);
        [VertexEdgeManifold1 NumberEdgeManifold1 NumberEdge1] = GetEdgeManifold(n,NumberVertex,LabelVertexManifold); 
        ne = 0;
        VertexEdgeManifold = [];
        for nk = 1:NumberEdge1
           if NumberEdgeManifold1(nk) > 2
               [Edge1 Edge2] = SolveAmbiguity(n,Grid,Pert,First,Delta,Func,LabelVertexManifold,VertexEdgeManifold1(nk,:)); 
               ne = ne+1; 
               VertexEdgeManifold = [VertexEdgeManifold; VertexEdgeManifold1(nk,Edge1(1)) VertexEdgeManifold1(nk,Edge1(2))];
               ne = ne+1; 
               VertexEdgeManifold = [VertexEdgeManifold; VertexEdgeManifold1(nk,Edge2(1)) VertexEdgeManifold1(nk,Edge2(2))];
           else
               ne = ne+1; 
               VertexEdgeManifold = [VertexEdgeManifold; VertexEdgeManifold1(nk,1:2)];
           end
        end
        NumberEdge = ne;
        if NumberVertex > 0
           fprintf(file,'%3d ',g);
           fprintf(file,'%3d ',Grid);
           fprintf(file,'\n');
           [nc nSkel Skel] = SkeletonHyperCube(file,n,NumberVertex,NumberEdge,VertexManifold,LabelVertexManifold,LabelVertex,VertexEdgeManifold);
           fprintf(file,'\n');
        end   
    
        for nck = 1:nc
            for i = 1:nSkel{nck}(n-k)
                [value g1] = NewHyperCube(n,Skel{nck}{n-k}(i),Grid,Division,Base);
                if value == 1
                    [value] = InList(NHcubeProc, HcubeProc, g1);
                    if value == 0
                        [NHcubeTrans HcubeTrans] = InsertList(NHcubeTrans,HcubeTrans,g1);
                    end
                end
            end
        end
   end
   fprintf(file,'-1\n');
   
   fclose(file);
   return
end


function [value] = InList(NList, List, g)
    value = 0;
    if NList > 0
        for i = 1:NList
            if List(i,1) >= g
                if List(i,1) == g
                    value = 1;
                    return
                else
                    return
                end
            end
        end
    end
    return
end


function [g NList List] = GetAndRemoveList(NList, g, List)
    if NList == 0
        g  = -1;
        return
    end
    for i = 1:NList
        if List(i,1) == g
            List  = [List(1:i-1,:); List(i+1:NList,:)];
            NList = NList - 1;
            return
        end
    end
    g     = List(1,1);
    List  = List(2:NList,:);
    NList = NList - 1;
    return
end


function [NList List] = InsertList(NList, List, g)
    if NList > 0
        for i = 1:NList
            if List(i,1) >= g
                if List(i,1) == g
                    return
                else
                    List  = [List(1:i-1,:); g; List(i:NList,:)];
                    NList = NList + 1;
                    return
                end
            end
        end
    end
    List  = [List; g];
    NList = NList + 1;
    return
end



function [value g] = NewHyperCube(n, i, Grid, Division, Base)
    if i >= n
        i = i - n;
        if Grid(i+1) < Division(i+1)-1
            G       = Grid;
            G(i+1)  = G(i+1) + 1;
            g       = Get_Label_Prod(n,Base,G);
            value   = 1;
            return
        else
            g       = 0;
            value   = 0;
            return
        end
    else
        if Grid(i+1) > 0
            G       = Grid;
            G(i+1)  = G(i+1) - 1;
            g       = Get_Label_Prod(n,Base,G);
            value   = 1;
            return
        else
            g       = 0;
            value   = 0;
            return
        end
    end
    return
end
          



function [value g s] = Pivoting(n, i, Grid, P, Division, Base, BaseP)
    if i == 0
        if Grid(P(1)) < Division(P(1))-1
            G       = Grid;
            G(P(1)) = G(P(1)) + 1;
            g       = Get_Label_Prod(n,Base,G);
            Q       = [P(2:n) P(1)];
            s       = Get_Label_Perm(n,BaseP,Q);
            value   = 1;
            return
        else
            g       = 0;
            s       = 0;
            value   = 0;
            return
        end
    end
    if i == n
        if Grid(P(n)) > 0
            G       = Grid;
            G(P(n)) = G(P(n)) - 1;
            g       = Get_Label_Prod(n,Base,G);
            Q       = [P(n) P(1:n-1)];
            s       = Get_Label_Perm(n,BaseP,Q);
            value   = 1;
            return
        else
            g       = 0;
            s       = 0;
            value   = 0;
            return
        end
    end
    g     = Get_Label_Prod(n,Base,Grid);
    Q     = [P(1:i-1) P(i+1) P(i) P(i+2:n)];
    s     = Get_Label_Perm(n,BaseP,Q);
    value = 1;
    return
end
          

function [g value] = GetFirstHypercube(n, k, First, Last, Division, FirstPoint, Func)
   tol  = 0.000001;
   [F]  = Func(n,FirstPoint);
   if norm(F(1:k),2) > tol
       g     = -1;
       s     = -1;
       value = -1;
       fprintf('Ponto nao pertence a variedade\n');
       return
   end
   for i = 1:n
        if (FirstPoint(i) < First(i)) || (FirstPoint(i) > Last(i))
            g     = -1;
            s     = -1;
            value = -2;
            fprintf('Ponto fora do dominio\n');
            return
        end
   end 
   Delta = (Last-First)./Division;
   Pert  = random('Normal',0,1,2^n,n)*tol;
   Base  = InitGrid(n, Division);
   for i = 1:n
        DivisionP(i) = i;
   end
   [BaseP] = InitGrid(n,DivisionP);    
   Grid    = fix((FirstPoint-First)./Delta);
   g    = Get_Label_Prod(n,Base,Grid);
   s    = 0;
   gold = -1; 
   while g > 0
        fprintf('Simplexo: g = %d  s = %d\n',g,s);
        if g ~= gold
            [Grid]       = GridCoords(n,g,Base);
            [Vert FVert] = GenVert(n, Grid, Pert, First, Delta, Func);
        end
        [P]     = Get_Perm(n,BaseP,s);
        [Simp]  = GenLabelSimplex(n,P); 
        [trans] = GetBestSimplex(n,Simp,Vert,FirstPoint);
        if trans < 0
            value = 1;
            return
        end 
        gold = g;
        [value g1 s1] = Pivoting(n,trans,Grid,P,Division,Base,BaseP);
        if value == 0
            g  = -1;
        else
            g  = g1;
            s  = s1;
        end 
   end
   g     = -1;
   s     = -1;
   value = -2;
   fprintf('Ponto fora do dominio\n');
   return
end


function [trans] = GetBestSimplex(n,Simp,Vert,FirstPoint)        
   for i = 1:n+1
      A(1,i) = 1;
      A(2:n+1,i) = Vert(Simp(i)+1,:);
   end
   b = [1; FirstPoint(1:n)'];
   lamb = A\b;
   trans = 0;
   if lamb >= 0.0    
       trans = -1;
   else
       min  = lamb(1);
       imin = 1;
       for i = 1:n+1
           if lamb(i) < min
               min  = lamb(i);
               imin = i;
           end
       end
       trans = imin-1;
   end   
   return
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


function [V FV] = GenVertHyperCube(n,Grid,Pert,First,Delta,Label,Func)     
   [Coords]    = HyperCubeCoords(n,Label);
   [CoordPert] = HyperCubePert(n,Grid,Coords,Pert);
   [V]         = HyperCube(n,First,Delta,Grid,Coords,CoordPert);
   [FV]        = Func(n,V); 
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










function [i] = GetPivotSimplex(n, Simp, Face)
    for i = 1:n
        vertex = Simp(i);
        [k]    = Include(1,vertex,n-1,Face);
        if k == 0
            return
        end
    end
    i = 0;
    return
end

function [f] = Get_Prod(n, Base, labelp)
    [f] = GridCoords(n,labelp,Base);
    return
end

function [p] = Get_Perm(n, Base, labelp)
    labelp = Base(n+1)-1-labelp;
    [f]    = GridCoords(n,labelp,Base);
    f      = f + 1;
    [p]    = Map_Perm(n,f);
    return
end

function [label] = Get_Label_Perm(n, Base, p)
    [f]     = Map_Inv_Perm(n,p);
    [label] = Get_Label_Prod(n,Base,f);
    label   = Base(n+1)-1 - label;
    return
end

function [label] = Get_Label_Prod(n, Base, f)
    label = 0;
    for i = 1:n
        label = label + f(i)*Base(i);
    end
    return
end

function [f] = Map_Inv_Perm(n,p)
    q = p;
    for i = n:-1:2
        [j]  = Find(i,i,q);
        f(i) = j;
        q    = [q(1:j-1) q(j+1:i)];
    end
    f(1) = 1;
    f    = f - 1;
    return
end

function [j] = Find(n,i,x)
    j = 0;
    for k = 1:n
        if x(k) == i
            j = k;
            return
        end
    end
    return
end


function [nc nSkel Skel] = SkeletonHyperCube(file, n, nv, ne, Vertex, FaceVertex, LabelVertex, Edges)
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
        [nSkel{k} Skel{k} AdjSkel{k}] = SkeletonComp(n, nVertex, FaceVertexComp);
        for j = 2:n
            fprintf(file,'%3d\n',nSkel{k}(j));
            for i = 1:nSkel{k}(j)
                fprintf(file,'%3d ',AdjSkel{k}{j}{i}(:));
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