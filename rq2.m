function [R, Q] = rq2(M)
    [Q,R] = qr(flipud(M)');
    R = flipud(R');
    R = fliplr(R);

    Q = Q';   
    Q = flipud(Q); 
end