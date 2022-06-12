
using Printf


function write_n_variables(file_name, arr1, num, dim)
  szz = size(arr1);

  nonZero_cnt = 0

  open(file_name, "w") do fileID
    if dim == 3
      for ii = 1:szz[1]
        for jj = 1:szz[2]
          for kk = 1:szz[3]
            val = arr1[ii,jj,kk]
            if val > 0
              nonZero_cnt = nonZero_cnt + 1
            
              write(fileID, string(val))
              write(fileID, "\n")

              if nonZero_cnt == num
                return
              end
            end
          end
        end
      end
    end

    if dim == 2
      
      for ii = 1:szz[1]
        for jj = 1:szz[2]
            val = arr1[ii,jj]
            if val > 0
              nonZero_cnt = nonZero_cnt + 1
            
              write(fileID, string(val))
              write(fileID, "\n")

              if nonZero_cnt == num
                return
              end
          end
        end
      end
    end


    if dim == 1
      for ii = 1:szz[1]
        val = arr1[ii]
        if val != 0
          nonZero_cnt = nonZero_cnt + 1
          write(fileID, string(val))
          write(fileID, "\n")
          if nonZero_cnt == num
            return
          end
        end
      end
    end
  end

  @printf("FILE %s : Non zero count = %d\n", file_name, nonZero_cnt);
  
end



function write_variable(file_name, arr1, dim)
  szz = size(arr1);

  nonZero_cnt = 0

  open(file_name, "w") do fileID
    if dim == 3
      for ii = 1:szz[1]
        for jj = 1:szz[2]
          for kk = 1:szz[3]
            val = arr1[ii,jj,kk]
            if val != 0
              nonZero_cnt = nonZero_cnt + 1
            end
            write(fileID, string(val))
            write(fileID, "\n")
          end
        end
      end
    end

    if dim == 1
      for ii = 1:szz[1]
        val = arr1[ii]
        if val != 0
          nonZero_cnt = nonZero_cnt + 1
        end
        write(fileID, string(val))
        write(fileID, "\n")
      end
    end
    write(fileID, " Non zero count = ")
    write(fileID, string(nonZero_cnt))
    write(fileID, "\n")
  end

  @printf("FILE %s : Non zero count = %d\n", file_name, nonZero_cnt);
  
end



function write_bool_variable(file_name, arr1, dim)
  szz = size(arr1);
  open(file_name, "w") do fileID

    if dim == 3
      for ii = 1:szz[1]
        for jj = 1:szz[2]
          for kk = 1:szz[3]
            if arr1[ii,jj,kk] == true
               write(fileID, string(1.0))
            else
              write(fileID, string(0.0))
            end
            write(fileID, "\n")
          end
        end
      end
    end
  end
end