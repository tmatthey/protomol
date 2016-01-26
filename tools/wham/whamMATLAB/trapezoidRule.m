function result = trapezoidRule(x,y)
    result = 0;
    
    for i=2:length(x)
        result = result + 0.5*(x(i) - x(i-1))*(y(i) + y(i-1));
    end