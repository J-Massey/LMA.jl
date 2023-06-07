function EUL_to_LAG(ni, nj, nk, nr, grdnw, varnw, delT, fn, Nt, omL, varsn, volsn)
    omLL = copy(omL)
    varsnL = copy(varsn)
    volsnL = copy(volsn)
    
    # Initialize omL
    omL[:, :, :, 1, 1:3] .= grdnw[:, :, :, 1:3]
    omLL[:, :, :, 1:3] .= omL[:, :, :, 1, 1:3]
    varsn[:, :, :, 1, 1:nr] .= varnw[:, :, :, 1, 1:nr]
    varsnL[:, :, :, 1:nr] .= varsn[:, :, :, 1, 1:nr]
    volsn[:, :, :, 1, 1] .= volsn[:, :, :, 1, 2]
    volsnL[:, :, :] .= volsn[:, :, :, 1, 1]
    
    i = 1
    j = 1
    k = 1
    
    for iNt = 2:Nt
        # Calculate omL surface position
        omLL[:, :, :, 1:3] .+= delT .* varsnL[:, :, :, 1:3]
        
        # Calculate velocity on omL surface
        for ii in 1:ni, jj in 1:nj, kk in 1:nk
            for i = 1:ni-1
                if omLL[ii, jj, kk, 1] ≥ grdnw[i, j, k, 1] && omLL[ii, jj, kk, 1] < grdnw[i+1, j, k, 1]
                    break
                end
            end
            for j = 1:nj-1
                if omLL[ii, jj, kk, 2] ≥ grdnw[i, j, k, 2] && omLL[ii, jj, kk, 2] < grdnw[i, j+1, k, 2]
                    break
                end
            end
            for k = 1:nk-1
                if omLL[ii, jj, kk, 3] ≥ grdnw[i, j, k, 3] && omLL[ii, jj, kk, 3] < grdnw[i, j, k+1, 3]
                    break
                end
            end
            if fn == 1
                varsnL[ii, jj, kk, 1:nr] = varnw[i, j, k, 1, 1:nr]
                volsnL[ii, jj, kk] = volsn[i, j, k, 1, 2]
            else
                varsnL[ii, jj, kk, 1:nr] = varnw[i, j, k, iNt, 1:nr]
                volsnL[ii, jj, kk] = volsn[i, j, k, iNt, 2]
            end
        end
        
        # Store omL and varsn
        omL[:, :, :, iNt, 1:3] .= omLL[:, :, :, 1:3]
        varsn[:, :, :, iNt, 1:nr] .= varsnL[:, :, :, 1:nr]
        volsn[:, :, :, iNt, 1] .= volsnL[:, :, :] ./ sum(volsnL[:, :, :])
        
        println(iNt, "/", Nt)
    end
end

function test_EUL_to_LAG()
    # Test data
    ni = 10
    nj = 10
    nk = 10
    nr = 5
    Nt = 3
    fn = 2
    delT = 0.1
    
    grdnw = rand(ni, nj, nk, 3)
    varnw = rand(ni, nj, nk, fn, nr)
    omL = zeros(ni, nj, nk, Nt, 3)
    varsn = zeros(ni, nj, nk, Nt, nr)
    volsn = ones(ni, nj, nk, Nt, 2)
    
    # Run the function
    EUL_to_LAG(ni, nj, nk, nr, grdnw, varnw, delT, fn, Nt, omL, varsn, volsn)
    
    # Assertions or other validation
    # ...
    println("Test passed!")
end

# Run the test
test_EUL_to_LAG()
