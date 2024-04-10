
function  [sig,pval] = perm_1d_delta(data1,data2,iter)

 

    size1 = size(data1,1);

    size2 = size(data2,1);

   

    null_delta = zeros(iter,1);

   

    data_combo = vertcat(data1,data2);

    grp_combo(1:size1,1) = 1;

    grp_combo(size1+1:size1+size2,1) = 2;

    orig_delta = nanmean(data_combo(grp_combo==1))-nanmean(data_combo(grp_combo==2));


    for x = 1:iter

        rand_vec = rand(size1+size2,1);

        [~,sort_rand] = sort(rand_vec);

        grp_rand = grp_combo(sort_rand);

        null_delta(x,1) = mean(data_combo(grp_rand==1))-mean(data_combo(grp_rand==2));

    end

   

    pos_thr = prctile(null_delta,95);

    sig = double(orig_delta > pos_thr);

    pval = sum(orig_delta>pos_thr)./iter;

   

end

 