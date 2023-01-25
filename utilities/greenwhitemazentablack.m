
function A = greenwhitemazentablack(ncolor)

    A = zeros(ncolor*2+2,3);
    A(1,1:3) = 0;
    A(2:ncolor+2,1) = (0:ncolor)/ncolor;
    A(2:ncolor+2,2) = 1;
    A(2:ncolor+2,3) = (0:ncolor)/ncolor;
    A(ncolor+2:2*ncolor+2,1) = 1;
    A(ncolor+2:2*ncolor+2,2) = (ncolor:-1:0)/ncolor;
    A(ncolor+2:2*ncolor+2,3) = 1;

end