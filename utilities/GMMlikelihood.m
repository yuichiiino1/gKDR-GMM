function L = GMMlikelihood(source_train, target_train, gm, B)

%%%%%%%%%%%  perform likelihood estimation  %%%%%%%%%%%%%%

assert(size(source_train,1) ~= 0);
assert(size(source_train,1) == size(target_train,1))
assert(size(target_train,2) == 1)

x_test = source_train*B;

L = zeros(size(x_test,1),1);
for i = 1:size(x_test,1)
    xi = x_test(i,:);
    gm_xi = conditional_GMM(gm, xi);
    L(i) = pdf(gm_xi, target_train(i));
end

end
