%[vert,tri] = icosphere(2);
[tri,vert] = stanford_bunny();
vert = vert./max(vert(:));

colormap(flip(winter))
h = trisurf(tri, vert(:,1), vert(:,2), vert(:,3)); %p(:,k)
daspect([1 1 1])
camlight
drawnow
%%

[L,M] = get_FEM_Lap(vert,tri);
%%
subplot(2,3,1)
spy(L)
title('S')

p = symrcm(L);
subplot(2,3,2)
spy(L(p,p))
title('S(p,p) R. C-M') %Cuthill-McKee

subplot(2,3,3)
spy(chol(L(p,p)))
title('chol(S(p,p)) R. C-M')

subplot(2,3,4)
spy(L)
title('S')

p = symamd(L);
subplot(2,3,5)
spy(L(p,p))
title('S(p,p) Sym. AMD') % approximate minimum degree permutation

subplot(2,3,6)
spy(chol(L(p,p)))
title('S(p,p) Sym. AMD')

function [L,M] = get_FEM_Lap(vert,tri)
    ntris = size(tri,1);
    nvert = size(vert,1);
    
    row = [];
    col = [];
    L_val = [];
    M_val = [];
    
    for k = 1:ntris
        t = tri(k,:);
        r1  = vert(t(1),:); %vertex i
        r2  = vert(t(2),:); %vertex j
        r3  = vert(t(3),:); %vertex l
        
        E3 = r1-r2; % edge 1
        E2 = r3-r1; % edge 2
        E1 = r2-r3; % edge 3
        
        A_tri = 0.5*norm(cross(E1,E2)); % area of tri
        
        % set local stiffnees
        A_k = (1.0/(4.0*A_tri)) * [[dot(E1,E1) dot(E1,E2) dot(E1,E3)];
                                   [dot(E2,E1) dot(E2,E2) dot(E2,E3)];
                                   [dot(E3,E1) dot(E3,E2) dot(E3,E3)]];
        
        % set local mass
        M_k = (1.0/12.0)*A_tri * [[2 1 1];
                                  [1 2 1];
                                  [1 1 2]];
                              
        % set entries of global mass and stiffness matricees
        for rr = 1:3
            for cc = 1:3
                row(end+1) = t(rr);
                col(end+1) = t(cc);
                L_val(end+1) = A_k(rr,cc);
                M_val(end+1) = M_k(rr,cc);
            end
        end
                 
    end
    
    L = sparse(row,col,L_val,nvert,nvert);
    M = sparse(row,col,M_val,nvert,nvert);
    
end