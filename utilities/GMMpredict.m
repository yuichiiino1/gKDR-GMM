function y = GMMpredict(source_test, gm, B)

%%%%%%%%%%%  perform prediction  %%%%%%%%%%%%%%

assert(size(source_test,1)~=0);

max_value = 3;

x_test = source_test*B;
u = x_test;

y = zeros(size(x_test,1),1);
for i = 1:size(x_test,1)
    xi = x_test(i,:);
    gm_xi = conditional_GMM(gm, xi);
    %disp(size(gm_xi.mu));
    %disp(size(gm_xi.Sigma));
    for count = 1:100
        Y = random(gm_xi);
        if Y>= -max_value && Y<=max_value
            break
        else
            Y = 0;
        end   
    end
    y(i) = Y;
end


end

