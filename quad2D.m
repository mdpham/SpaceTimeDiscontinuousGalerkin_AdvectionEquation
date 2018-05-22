function approx = quad2D(fun)
    s3 = sqrt(1/3);
    quadpts = [[-s3,-s3];
        [s3,-s3];
        [-s3,s3];
        [s3,s3]];
    quadwts = 1;
    approx = 0;
    for i = 1:4
        xy = quadpts(i,:);
        x = xy(1); y = xy(2);
        temp = fun(x,y)*quadwts;
        approx = approx+temp;
    end
end