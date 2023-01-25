function y = GPpredict_det(source_test, gp, B)

%%%%%%%%%%%  perform prediction  %%%%%%%%%%%%%%

assert(size(source_test,1)~=0);

max_value = 3;

x_test = source_test*B;
%u = x_test;

y = zeros(size(x_test,1),1);
for i = 1:size(x_test,1)
    xi = x_test(i,:);
    [Y,~,~] = predict(gp,xi); % [ypred,ysd,yint]
    
    if Y < -max_value || Y > max_value
        Y = 0;
    end
    
    y(i) = Y;
end

end

