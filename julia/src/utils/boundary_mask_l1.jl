

function boundary_mask_l1(B, G, sz, l1)

    nx = sz[0];
    ny = sz[1];
    nz = sz[2];
    nxny = nx*ny;
    nxnynz = nx*ny*nz;

    NX = nx-1;
    NY = nx*(ny-1);
    NZ = nxny*(nz-1)

    return

    #uint8_t *GL = (uint8_t *)calloc(nx*ny*nz, sizeof(*G));
    GL = zeros(nx*ny*nz)

    # pragma omp parallel for private(i,j,k,l) schedule(static) \
    #    if(nxny*nz > 32*32*32)
    for k = 1:nxnynz 
        for j = 1:nxny
            l = j + k
            for i=1:nx
                if ((i == 0) || (j == 0) || (k == 0)  || (i == NX) || (j == NY) || (k == NZ)) 
                    B[l] = GL[l] #; // = G[l];
                else
                    #GL[l] = G[l] && !(G[l-1] && G[l+1] && G[l-nx] && G[l+nx] && G[l-nxny] && G[l+nxny]
                    B[l] = GL[l] 
                end
                i = i + 1
                l = l + 1
            end
            j += nx
        end
        
        k += nxny
    end

    return
end
