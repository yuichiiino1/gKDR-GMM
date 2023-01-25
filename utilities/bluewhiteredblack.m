
function A = bluewhiteredblack(ncolor)

    A = zeros(ncolor*2+2,3);
    A(2:ncolor+2,1) = (0:ncolor)/ncolor;
    A(2:ncolor+2,2) = (0:ncolor)/ncolor;
    A(2:ncolor+2,3) = 1;
    A(ncolor+2:2*ncolor+2,1) = 1;
    A(ncolor+2:2*ncolor+2,2) = (ncolor:-1:0)/ncolor;
    A(ncolor+2:2*ncolor+2,3) = (ncolor:-1:0)/ncolor;

end