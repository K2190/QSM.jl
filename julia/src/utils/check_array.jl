
function check_array(x)

    sz_m0 = size(x);
    count_m0 = 0;
    for ii=1:sz_m0[1]
        for jj=1:sz_m0[2]
            for kk=1:sz_m0[3]
                val = x[ii, jj, kk]
                if val > 0
                    count_m0 = count_m0 + 1;
                end
            end
        end
    end
    println("Array is greater than zero for ")
    println(count_m0)
    
    count_m0 = 0;
    for ii=1:sz_m0[1]
        for jj=1:sz_m0[2]
            for kk=1:sz_m0[3]
                val = x[ii, jj, kk]
                if val < 0
                    count_m0 = count_m0 + 1;
                end
            end
        end
    end
    println("Array is less than zero for ");
    println(count_m0)
    
    count_m0 = 0;
    for ii=1:sz_m0[1]
        for jj=1:sz_m0[2]
            for kk=1:sz_m0[3]
                val = x[ii, jj, kk]
                if val == 0
                    count_m0 = count_m0 + 1;
                end
            end
        end
    end
    println("Array is equal to zero for ")
    println(count_m0)
    
end    