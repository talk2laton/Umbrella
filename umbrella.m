
close all;
hg1 = NewFigure;

[Xc, Yc, Zc] = cylinder([0.02,0.04,0.04], 20);
Zc(1,:) = 5.5; Zc(2,:) = 5; Zc(3,:) = 1;
surf(0.25 + Xc, Yc, Zc, EdgeAlpha = 0, FaceColor = 'b',Parent = hg1); 
view(3); axis equal; hold on;
axis([-4,4,-4,4,-1,6]);

[Xs, Ys, Zs] = cylinder([0.04,0.07,0.07,0.05,0.05], 20);
Zs([1,2],:) = 0.1; Zs([3,4],:) = 0.09; Zs(5,:) = 0.0;
Xs = 0.25 + [Xs; flipud(Xs)]; Ys = [Ys; flipud(Ys)]; 
Zs = 3.9 + [Zs; -flipud(Zs)];
surf(Xs, Ys, Zs, FaceColor = 'k', Parent = hg1);
axis([-4,4,-4,4,-1,6]);
[Xh, Yh, Zh] = cylinder([0.04,0.06,0.06], 20);
Zh(1,:) = 1; Zh(2,:) = 0.9; Zh(3,:) = 0; Xh = 0.25 + Xh;
x = Xh(end,:); y = Yh(end,:); z = Zh(end,:);
for i = 1:12
    [x, y, z] = rotate(x, y, z, [0,1,0], pi/12);
    Xh = [Xh; x]; Yh = [Yh; y]; Zh = [Zh; z];
end
surf(Xh, Yh, Zh, EdgeAlpha = 0.2, FaceColor = 'g', ...
    FaceAlpha = 0.8, Parent = hg1);

t = linspace(0,pi/4,51)';
u = pi/12*linspace(-1, 1, 51); 
R = 5+2*t*u.^2;
X = R.*(sin(t)*cos(u)); 
Y = R.*(sin(t)*sin(u)); 
Z = R.*(cos(t)*ones(size(u))); 
x = 0.07*cos(pi/12*(1:2:23)); 
y = 0.07*sin(pi/12*(1:2:23));
for i = 1:12
    xyzl = [X(:,end), Y(:,end), Z(:,end)];
    xyzl = [xyzl; 2*xyzl(end,:)-xyzl(end-1,:)];
    xl = xyzl(:,1); yl = xyzl(:,2); zl = xyzl(:,3);
    surf(0.25 + X, Y, Z, 'EdgeAlpha',0, FaceColor = 'r', Parent = hg1); 
    plot3(0.25+[x(i), xl(25)], [y(i), yl(25)], [4, zl(25)], 'k', ...
        LineWidth = 1.5, Parent = hg1);
    plot3(0.25 + xl, yl, zl, 'k', LineWidth = 1.5, Parent = hg1);
    plot3(0.25 + X(end,:), Y(end,:), Z(end,:), 'k', Parent = hg1);
    [X, Y, Z] = rotate(X, Y, Z, [0,0,1], pi/6);
 end

axis([-4,4,-4,4,-1,6]); axis equal
saveas(gcf,'pic1.png')
set(hg1,'Matrix', makehgtform('yrotate',pi/3));
saveas(gcf,'pic2.png')
set(hg1,'Matrix', makehgtform('xrotate',-pi));
saveas(gcf,'pic3.png')
set(hg1,'Matrix', makehgtform('yrotate',-pi/3));
saveas(gcf,'pic4.png')
set(hg1,'Matrix', makehgtform('xrotate',-pi/3));
saveas(gcf,'pic5.png')
im1 = imread('pic1.png'); im2 = imread('pic2.png');
im3 = imread('pic3.png'); im4 = imread('pic4.png');
im5 = imread('pic5.png'); 
imwrite([im1, im2, im3, im4, im5], 'pic6.png');

for i = 1:800
    set(hg1,'Matrix', makehgtform('zrotate',i*pi/400));
    drawnow
end

function hg1 = NewFigure
	figure(Color = 'w', Position = [680 54 1021 924])
	ax = gca; hg1 = hgtransform('Parent',ax); 
	camlight('headlight', 'infinite '); lighting gouraud;
end