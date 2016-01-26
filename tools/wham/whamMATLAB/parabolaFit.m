function result = parabolaFit(x,y)
    result = [];
    for i = 1:(length(x) - 1)
        if i == 1
            % Special case
            pt1 = 1;
            pt2 = 2;
            pt3 = 3;
            first = x(1);
            last = x(2);
        else
            pt1 = i - 1;
            pt2 = i;
            pt3 = i + 1;
            first = x(i);
            last = x(i + 1);
        end
        
        % Parabola coefficients
        X = [x(pt1)*x(pt1) x(pt1) 1;
             x(pt2)*x(pt2) x(pt2) 1;
             x(pt3)*x(pt3) x(pt3) 1];
        Y = [y(pt1);
             y(pt2);
             y(pt3)];
        A = X \ Y;
        
        newx = transpose(first:0.01:last);
        newy = polyval(A,newx);
        new = [newx newy];
        result = [result; new];
    end