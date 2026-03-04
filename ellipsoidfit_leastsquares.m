function p = ellipsoidfit_leastsquares(data)

%%%%2D ellipsoid fitting with LLS
%    Ax^2 + By^2 + Cxy + Dx + Ey + J = 0

%%%%3D ellipsoid fitting with LLS
%    the 10 parameters describing the ellipsoid algebraically
%    Ax^2 + By^2 + Cz^2 + Dxy + Exz + Fyz + Gx + Hy + Iz + J = 0
%    where norm([A,B,C,D,E,F,G,H,I,J]) == 1

dimension=size(data,2);
if dimension==2   % 2D ellipse
    x = data(:,1);
    y = data(:,2);
    assert(numel(x) >= 5, ...
    'At least 5 points are required to fit a unique ellipse.');
    % normalize data by canceling mean to improve numerical accuracy
    mx = mean(x);
    my = mean(y);
    x = x - mx;
    y = y - my;
    % use singular value decomposition (unconstrained problem)
    D = [ x .* x,  y .* y,  x .* y,  x, y, ones(size(x)) ];
    [U,S,V] = svd(D,0);
    p = V(:,end);  % smallest singular value
    % unnormalize
    p(:) = ...
    [ p(1); p(2); p(3) ...
    ; p(4) - 2*p(1)*mx - p(3)*my ...
    ; p(5) - p(3)*mx - 2*p(2)*my ...
    ; p(6) + p(1)*mx*mx + p(2)*my*my +  p(3)*mx*my  - p(4)*mx - p(5)*my ];
    p(3)=p(3)/2;
    p(4)=p(4)/2;
    p(5)=p(5)/2;
elseif dimension==3           % 3D ellipsoid
    x = data(:,1);
    y = data(:,2);
    z = data(:,3);  
    assert(numel(x) >= 9, ...
    'At least 9 points are required to fit a unique ellipsoid.');
    % normalize data by canceling mean to improve numerical accuracy
    mx = mean(x);
    my = mean(y);
    mz = mean(z);
    x = x - mx;
    y = y - my;
    z = z - mz;
    % use singular value decomposition (unconstrained problem)
    D = [ x .* x, y .* y,  z .* z, x .* y, x .* z, y .* z, x, y, z,  ones(size(x)) ];
    [U,S,V] = svd(D,0);
    p = V(:,end);  % smallest singular value
    % unnormalize
    p(:) =  [ p(1) ; p(2); p(3)  ; p(4) ; p(5) ; p(6) ...
    ; p(7) - 2*p(1)*mx - p(4)*my - p(5)*mz ...
    ; p(8) - 2*p(2)*my - p(4)*mx - p(6)*mz ...
    ; p(9) - 2*p(3)*mz - p(5)*mx - p(6)*my ...
    ; p(10) + p(1)*mx*mx + p(2)*my*my + p(3)*mz*mz + p(4)*mx*my + p(5)*mx*mz + p(6)*my*mz - p(7)*mx - p(8)*my - p(9)*mz ...
    ];
    p(4)=p(4)/2;
    p(5)=p(5)/2;
    p(6)=p(6)/2;
    p(7)=p(7)/2;
    p(8)=p(8)/2;
    p(9)=p(9)/2;
else
    error('No LLS in high-dimensional');
end
% normalize
p = p / norm(p,2);
if p(1)<0
    p=-p;
end

