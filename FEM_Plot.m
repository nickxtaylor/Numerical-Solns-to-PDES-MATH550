load trimesh2d
close all
t = trife;
p = [xfe(:), yfe(:), 0*yfe(:)];
triplot(trife,xfe,yfe)
h = trisurf(t, p(:,1), p(:,2), p(:,3));
set(h,'linewidth',3,'facecolor',[1 1 1])
hold all

II = 667; %randi(length(p),1)

%cols = {'c','m','y'};
cols = {'r','g','b'};

for ijk = 1:3
    v = t(II,ijk);
    f = p;
    f(:,3) = -0.01;
    f(v,3) = 1.0;
    hh(ijk) = trisurf(t, f(:,1), f(:,2), f(:,3));
    set(hh(ijk),'linewidth',1,'edgecolor','k','facecolor',cols{ijk},'facealpha',0.2);
    hold all
end

xlim([0 60])
ylim([0 60])
view([-34.0875   59.6854])
%view([0   90])

%%
x = linspace(0,1,100);
[eta,xi] = meshgrid(x,x);

surf(xi,eta,max(1-eta-xi,0),'facecolor','r')
hold all
surf(xi,eta,xi.*(eta<(1-xi)),'facecolor','b')
xlabel('$$\xi$$')
ylabel('$$\eta$$')
camlight
