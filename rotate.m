function [X, Y, Z] = rotate(X, Y, Z, U, theta)
    M = Mxyz(U,theta);
    for i = 1:size(X,1)
        xyz= M*[X(i,:);Y(i,:);Z(i,:)];
        X(i,:) = xyz(1,:);
        Y(i,:) = xyz(2,:);
        Z(i,:) = xyz(3,:);
    end
end

function M = Mxyz(U,theta)
    c = cos(theta); s = sin(theta); ux = U(1); uy = U(2); uz = U(3);
    M = [ux*ux*(1 - c) +    c  uy*ux*(1 - c) - uz*s  uz*ux*(1 - c) + uy*s
         uy*ux*(1 - c) + uz*s  uy*uy*(1 - c) +    c  uz*uy*(1 - c) - ux*s
         uz*ux*(1 - c) - uy*s  uy*uz*(1 - c) + ux*s  uz*uz*(1 - c) +    c];  
end