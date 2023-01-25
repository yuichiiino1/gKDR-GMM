function y = GMMpredict_det(source_test, gm, B)

%%%%%%%%%%%  perform prediction  %%%%%%%%%%%%%%

assert(size(source_test,1)~=0);

max_value = 1;

x_test = source_test*B;
u = x_test;

y = zeros(size(x_test,1),1);
for i = 1:size(x_test,1)
    xi = x_test(i,:);
    gm_xi = conditional_GMM(gm, xi);
    
    
    
    %disp(size(gm_xi.mu));
    %disp(size(gm_xi.Sigma));
%     for count = 1:100
        %Y = random(gm_xi);%%%ここをgm_xiの期待値を計算に変わる。gm_xi.mu、gm_xi.Componentで、かけてわ。
        try 
            Y = gm_xi.ComponentProportion*gm_xi.mu;%gm_xi.muは3行１列、gm_xi.proportionは1行3列。
        catch
            save('error_det.mat','gm_xi')
        end
        
%         if Y>= -max_value && Y<=max_value
%             break
%         else
%             Y = 0;
%         end   


%         if (Y == 0)%ぴったし0になるのは、conditional_GMMでうまく作れなかったとき。このときpdfで当てはめて確率を出そうとすると、うまく行かない(分散が0だから。デルタ関数みたいな形の確率分布)
%             break
%         end
%         try 
%             if (pdf(gm_xi, Y) > 0.01)
%                 break
%             else
%                 Y = 0;
%             end
% 
%         catch
%             
%         end
        
    %end
    y(i) = Y;
    if (Y == 0)
        %disp('Cant make Y !!')
    end
end


end