function approx = quadrature2D(fun)
    quadpts = [[-0.57735,-0.57735];
        [0.57735,-0.57735];
        [-0.57735,0.57735];
        [0.57735,0.57735]];
    quadwts = 1;
    approx = 0;
    for i = 1:4
        xy = quadpts(i,:);
        x = xy(1); y = xy(2);
        temp = fun(x,y)*quadwts;
        approx = approx+temp;
    end
end