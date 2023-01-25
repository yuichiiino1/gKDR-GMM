
function A = bluewhitered(ncolor)

    A = zeros(ncolor*2+1,3);
    A(1:ncolor+1,1) = (0:ncolor)/ncolor;
    A(1:ncolor+1,2) = (0:ncolor)/ncolor;
    A(1:ncolor+1,3) = 1;
    A(ncolor+1:2*ncolor+1,1) = 1;
    A(ncolor+1:2*ncolor+1,2) = (ncolor:-1:0)/ncolor;
    A(ncolor+1:2*ncolor+1,3) = (ncolor:-1:0)/ncolor;

end