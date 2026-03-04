function model = Ours(data,lamda,stop_value)
    LLS_model = ellipsoidfit_leastsquares(data);
    %%%数据读取
    if size(data,2)==2      
        x=data(:,1);
        y=data(:,2);
        D = [ x .* x, y .* y,   2 * x .* y,  2 * x,2 * y, ones(size(x)) ];
        %定义模型方程和待求变量
        Q33=@(q) [ q(1), q(3); q(3), q(2)];
        q79=@(q) [q(4);q(5)];
        t=@(q) q79(q)'/Q33(q)*q79(q)-q(6); %
        r_vec=@(q) sqrt(t(q)./eig(Q33(q)));
        %%%%%投影距离1 d=s*r
        s=@(q) (sqrt(D*q/t(q)+1)-1);
        d_vec=@(q) s(q)*r_vec(q)'; %500*3矩阵
        d_projection=@(q) abs(sum(d_vec(q),2)/2);
        %%%% Sampson距离
        Fx=@(v) 2*v(1)*x+2*v(3)*y+2*v(4);
        Fy=@(v) 2*v(2)*y+2*v(3)*x+2*v(5);
        d_sampson =@(q) abs(D*q./sqrt(Fx(q).*Fx(q)+Fy(q).*Fy(q)));   %可正可负
    elseif size(data,2)==3
        x=data(:,1);
        y=data(:,2);
        z=data(:,3);
        D = [ x .* x, y .* y,  z .* z, 2 * x .* y,  2 * x .* z, 2 * y .* z,  2 * x,2 * y, 2 * z,ones(size(x)) ];
        %定义模型方程和待求变量
        Q33=@(q) [ q(1) q(4) q(5); ...
                   q(4) q(2) q(6); ...
                   q(5) q(6) q(3)];
        q79=@(q) [q(7);q(8);q(9)];
        t=@(q) q79(q)'/Q33(q)*q79(q)-q(10); %
        r_vec=@(q) sqrt(t(q)./eig(Q33(q)));
        %%%%%投影距离1 d=s*r
        s=@(q) (sqrt(D*q/t(q)+1)-1);
        d_vec=@(q) s(q)*r_vec(q)'; %500*3矩阵
        d_projection=@(q) abs(sum(d_vec(q),2)/3);
        %%%% Sampson距离
         Fx=@(q) 2*q(1)*x+2*q(4)*y+2*q(5)*z+2*q(7);
         Fy=@(q) 2*q(2)*y+2*q(4)*x+2*q(6)*z+2*q(8);
         Fz=@(q) 2*q(3)*z+2*q(5)*x+2*q(6)*y+2*q(9);
         d_sampson =@(q) abs(D*q./sqrt(Fx(q).*Fx(q)+Fy(q).*Fy(q)+Fz(q).*Fz(q)));   %可正可负
    else
         error('Dimension is wrong');
    end
    %%%%距离联合
    d_combinition=@(q) lamda*d_projection(q)+(1-lamda)*d_sampson(q);
    cost=@(q) d_combinition(q)+max(eig(Q33(q))<(10e-7))*100;
    [q, s, cnt] = LMFnlsq(cost,LLS_model,'Display',1,'FunTol',stop_value,'XTol',stop_value);
    q = q / norm(q,2);
    if q(1)<0
        q=-q;
    end
    model=q;
end
