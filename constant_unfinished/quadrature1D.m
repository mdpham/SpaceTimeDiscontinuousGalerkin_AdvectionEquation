% FIVE POINT QUADRATURE
function approx = quadrature1D(fun)
    quadpts = [0,...
        sqrt(5-2*sqrt(10/7))/3, -sqrt(5-2*sqrt(10/7))/3,...
        sqrt(5+2*sqrt(10/7))/3, -sqrt(5+2*sqrt(10/7))/3];
    quadwts = [128/225,...
        (322+13*sqrt(70))/900, (322+13*sqrt(70))/900,...
        (322-13*sqrt(70))/900, (322-13*sqrt(70))/900];
    approx = 0;
    for i = 1:5
        temp = fun(quadpts(i))*quadwts(i);
        approx = approx+temp;
    end
end