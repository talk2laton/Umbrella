function Umbrella = MakeUmbrella(ax)
if(nargin == 0)
    figure(Color = 'w');
    ax = gca; 
end
camlight(ax, 'headlight', 'infinite ');
hg1 = hgtransform('Parent',ax); 
[Xc, Yc, Zc] = cylinder([0.02,0.04,0.04], 20);
Zc(1,:) = 5.5; Zc(2,:) = 5; Zc(3,:) = 1;
surf(0.25 + Xc, Yc, Zc, EdgeAlpha = 0, FaceColor = 'b',Parent = hg1); 
view(3); axis equal; hold on;
fill3(0.25 + 7/4*Xc(2,:), 7/4*Yc(2,:), Zc(2,:),'g');

[Xs, Ys, Zs] = cylinder([0.04,0.07,0.07,0.05,0.05], 20);
Zs([1,2],:) = 0.1; Zs([3,4],:) = 0.09; Zs(5,:) = 0.0;
Xs = 0.25 + [Xs; flipud(Xs)]; Ys = [Ys; flipud(Ys)]; 
Zs = 3.9 + [Zs; -flipud(Zs)];
Slive = surf(Xs, Ys, Zs, FaceColor = 'k', Parent = hg1);
axis([-4,4,-4,4,-1,6]);

[Xh, Yh, Zh] = cylinder([0.04,0.06,0.06], 20);
Zh(1,:) = 1; Zh(2,:) = 0.9; Zh(3,:) = 0; Xh = 0.25 + Xh;
Handle = surf(Xh, Yh, Zh, EdgeAlpha = 0, FaceColor = 'g', Parent = hg1);
x = Xh(end,:); y = Yh(end,:); z = Zh(end,:);
for i = 1:12
    [x, y, z] = rotate(x, y, z, [0,1,0], pi/12);
    Xh = [Xh; x]; Yh = [Yh; y]; Zh = [Zh; z];
	Handle.XData = Xh;  Handle.YData = Yh;  Handle.ZData = Zh; 
end

t0 = fzero(@(t)(5 + 5*t)*sin(t)-0.07, 0);
t = linspace(t0,0.8,21)'; u = pi/12*linspace(-1, 1, 21);  R = 5 + 5*t*u.^2;
X = R.*(sin(t)*cos(u)); Y = R.*(sin(t)*sin(u)); Z = 5*(cos(t)*ones(size(u)));
x = 0.07*cos(pi/12*(1:2:23)); y = 0.07*sin(pi/12*(1:2:23)); 
Xv = X; Yv = Y; Zv = Z;
xyzl = [X(:,end), Y(:,end), Z(:,end)];
xyzl = [xyzl; 1.5*xyzl(end,:)-0.5*xyzl(end-1,:)];
xl = xyzl(:,1); yl = xyzl(:,2); zl = xyzl(:,3);
Xl = xl; Yl = yl; Zl = zl;
Xp = [x(1), xl(10)]; Yp = [y(1), yl(10)]; Zp = [4, zl(10)];
for i = 2:12
    [X, Y, Z] = rotate(X, Y, Z, [0,0,1], pi/6);
    xyzl = [X(:,end), Y(:,end), Z(:,end)];
    xyzl = [xyzl; 1.5*xyzl(end,:)-0.5*xyzl(end-1,:)];
    xl = xyzl(:,1); yl = xyzl(:,2); zl = xyzl(:,3);
    Xv = [Xv,X]; Yv = [Yv,Y]; Zv = [Zv,Z]; 
    Xl = [Xl; nan; xl]; Yl = [Yl; nan; yl]; Zl = [Zl; nan; zl];
    Xp = [Xp, nan, x(i), xl(10)]; Yp = [Yp, nan, y(i), yl(10)]; 
    Zp = [Zp, nan, 4, zl(10)];
end
Cover = surf(0.25 + Xv, Yv, Zv, 'EdgeAlpha',0, FaceColor = 'y', Parent = hg1); 
Frame = plot3(0.25 + Xl, Yl, Zl, 'k', LineWidth = 1.1, Parent = hg1);
Brace = plot3(0.25 + Xp, Yp, Zp, 'k', LineWidth = 1.1, Parent = hg1);
Edge = plot3(0.25 + Xv(end,:), Yv(end,:), Zv(end,:), 'k', Parent = hg1); 
axis([-4,4,-4,4,-1,6]); axis equal;

dxyzl = diff(xyzl);
L = cumsum([0; vecnorm(dxyzl, 2, 2)]);
dxyzv = diff([X(end,:);Y(end,:);Z(end,:)]');
V = cumsum([0; vecnorm(dxyzv, 2, 2)]);

CoverXOpen = Cover.XData;
CoverYOpen = Cover.YData;
CoverZOpen = Cover.ZData;

FrameXOpen = Frame.XData;
FrameYOpen = Frame.YData;
FrameZOpen = Frame.ZData;

BraceXOpen = Brace.XData;
BraceYOpen = Brace.YData;
BraceZOpen = Brace.ZData;

EdgeXOpen = Edge.XData;
EdgeYOpen = Edge.YData;
EdgeZOpen = Edge.ZData;

SliveXOpen = Slive.XData;
SliveYOpen = Slive.YData;
SliveZOpen = Slive.ZData;

t = linspace(-1,0.5)*pi/4;
r = cos(2*t);
x1 = -r.*cos(t); y1 = r.*sin(t);
t = linspace(-1,0)*pi/4; r = cos(2*t);
x2 = r.*cos(t); y2 = -r.*sin(t);
x = [0.2*flip(x1), x2]; y = [0.2*flip(y1), y2];
f = 0.2-0.1*(x+1); y = f+y;
x = [x, flip(x)]; y = [y, flip(-y)];
x([200,300]) = []; y([200,300]) = [];
l = cumsum([0,hypot(diff(x), diff(y))]);
l = V(end)/l(end)*l; l(end) = V(end);
x = interp1(l, x, V)';  y = interp1(l, y, V)';
x = x - min(x); y = flip(y);

r1 = 0.07;
r2 = 0.2;
dy = diff(y([1, end]));
X = linspace(r1,r2,21)'*ones(size(u)) + linspace(0,0.2,21)'*x;
Y = linspace(r1/dy,r2/dy,21)'*y/5;
Z = linspace(5, 0.9, 21)'*ones(size(u));
% x = X(:,1);y = Y(:,1);z = Z(:,1);
% dx = diff(x); dy = diff(y); dz = diff(z);
% lref = cumsum([0;sqrt(dx.^2 + dy.^2 + dz.^2)]);
% for n = 2:50
%     x = X(:,n);y = Y(:,n);z = Z(:,n);
%     dx = diff(x); dy = diff(y); dz = diff(z);
%     l = cumsum([0;sqrt(dx.^2 + dy.^2 + dz.^2)]);
%     X(:,n) = interp1(l, x, lref);  
%     Y(:,n) = interp1(l, y, lref);
%     Z(:,n) = interp1(l, z, lref);
% end
[Xf, Yf, Zf] = rotate(X, Y, Z, [0,0,1], pi/6);
[Xb, Yb, Zb] = rotate(X, Y, Z, [0,0,1], -pi/6);
X(:,end) = 0.5*(X(:,end) + Xf(:,1));
Y(:,end) = 0.5*(Y(:,end) + Yf(:,1));
Z(:,end) = 0.5*(Z(:,end) + Zf(:,1));

X(:,1) = 0.5*(X(:,1) + Xb(:,end));
Y(:,1) = 0.5*(Y(:,1) + Yb(:,end));
Z(:,1) = 0.5*(Z(:,1) + Zb(:,end));

x = 0.07*cos(pi/6*(0.5:12)); y = 0.07*sin(pi/6*(0.5:12)); 
Xv = X; Yv = Y; Zv = Z;
xyzl = [X(:,end), Y(:,end), Z(:,end)];
xyzl = [xyzl; 1.2*xyzl(end,:)-0.2*xyzl(end-1,:)];
xl = xyzl(:,1); yl = xyzl(:,2); zl = xyzl(:,3);
Xl = xl; Yl = yl; Zl = zl;
Xp = [x(1), xl(10)]; Yp = [y(1), yl(10)]; Zp = [1.3, zl(10)];
for i = 2:12
    [X, Y, Z] = rotate(X, Y, Z, [0,0,1], pi/6);
    xyzl = [X(:,end), Y(:,end), Z(:,end)];
    xyzl = [xyzl; 1.2*xyzl(end,:)-0.2*xyzl(end-1,:)];
    xl = xyzl(:,1); yl = xyzl(:,2); zl = xyzl(:,3);
    Xv = [Xv,X]; Yv = [Yv,Y]; Zv = [Zv,Z]; 
    Xl = [Xl; nan; xl]; Yl = [Yl; nan; yl]; Zl = [Zl; nan; zl];
    Xp = [Xp, nan, x(i), xl(10)]; Yp = [Yp, nan, y(i), yl(10)]; 
    Zp = [Zp, nan, 1.3, zl(10)];
end

CoverXClse = 0.25 + Xv;
CoverYClse = Yv;
CoverZClse = Zv;

FrameXClse = 0.25 + Xl';
FrameYClse = Yl';
FrameZClse = Zl';

BraceXClse = 0.25 + Xp;
BraceYClse = Yp;
BraceZClse = Zp;

EdgeXClse = 0.25 + Xv(end,:);
EdgeYClse = Yv(end,:);
EdgeZClse = Zv(end,:);

SliveXClse = SliveXOpen;
SliveYClse = SliveYOpen;
SliveZClse = SliveZOpen - 2.7;

Umbrella.Open = @(t) open(1-t);
Umbrella.Cover = Cover;
Umbrella.Edge = Edge;
    function open(t)
        Cover.XData =  CoverXOpen + t*(CoverXClse - CoverXOpen);
        Cover.YData =  CoverYOpen + t*(CoverYClse - CoverYOpen);
        Cover.ZData =  CoverZOpen + t*(CoverZClse - CoverZOpen);
    
        Frame.XData =  FrameXOpen + t*(FrameXClse - FrameXOpen);
        Frame.YData =  FrameYOpen + t*(FrameYClse - FrameYOpen);
        Frame.ZData =  FrameZOpen + t*(FrameZClse - FrameZOpen);
    
        Slive.XData =  SliveXOpen + t*(SliveXClse - SliveXOpen);
        Slive.YData =  SliveYOpen + t*(SliveYClse - SliveYOpen);
        Slive.ZData =  SliveZOpen + t*(SliveZClse - SliveZOpen);
    
        Brace.XData =  BraceXOpen + t*(BraceXClse - BraceXOpen);
        Brace.YData =  BraceYOpen + t*(BraceYClse - BraceYOpen);
        Brace.ZData =  BraceZOpen + t*(BraceZClse - BraceZOpen);
    
        Edge.XData =  EdgeXOpen + t*(EdgeXClse - EdgeXOpen);
        Edge.YData =  EdgeYOpen + t*(EdgeYClse - EdgeYOpen);
        Edge.ZData =  EdgeZOpen + t*(EdgeZClse - EdgeZOpen);
        drawnow;
    end
end



