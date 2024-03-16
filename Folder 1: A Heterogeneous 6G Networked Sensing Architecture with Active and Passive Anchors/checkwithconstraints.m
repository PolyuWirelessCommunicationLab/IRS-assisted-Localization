function da_final=checkwithconstraints(K,da_crrct,da_crrct_row)
%%This function aims to find all the data association solutions that satisfy (26) (27)
    for k = 1:K
        for i = 1:da_crrct_row
            if da_crrct(i,1) == k
                r(k+1) = i;
            end
        end
    end
    %%find all the possbile combination 
    da = {};
    for k = 1:K
        da{k} = [];
    end
    k_new = 1;
    while k_new <= K
        if k_new == 1
             for i = r(1)+1:r(2)
                da{1}(i,:) = da_crrct(i,1:6);
             end
        end
        if k_new >= 2
            [da_cur_row,da_cur_col] = size(da{k_new-1});
            for j = r(k_new)+1:r(k_new+1)
                for i = 1:da_cur_row 
                    flag = 0;
                    da_pre = da{k_new-1}(i,:);
                    for ind = 1:da_cur_col
                        index = mod(ind,6);
                        if index ~=0 &&index ~= 5
                            if da_crrct(j,index) == da_pre(ind)
                                flag = 1;
                            end
                        end
                    end
                    if flag ~= 1
                        da{k_new} = [da{k_new};da_pre,da_crrct(j,1:6)];
                    end
                end
            end 
        end
        k_new = k_new+1;
    end  
    [da_row,da_col] = size(da{K});
    da_final = [];
    for e = 1:da_row
        da_final = [da_final;reshape(da{K}(e,:),[6,K])'];  
    end
end
