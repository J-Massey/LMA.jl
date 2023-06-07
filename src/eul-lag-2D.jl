using Test
"""
    EUL_to_LAG(ni, nj, nk, nr, grdnw, varnw, delT, fn, Nt, omL, varsn, volsn)

Translate the reference frame of a vorticity field from Eulerian to Lagrangian.
The input to the algorithm is a set of Eulerian flow fields and associated flow
parameters that include the grid, the solution state variables, the number of snapshots
(Nt ) and the time step (δt)

## Arguments
- `ni`: Number of grid points in the x-direction.
- `nj`: Number of grid points in the y-direction.
- `nk`: Number of grid points in the z-direction.
- `nr`: Number of variables.
- `grdnw`: Array of shape `(ni, nj, nk, 3)` representing the gradient of the field.
- `varnw`: Array of shape `(ni, nj, nk, fn, nr)` representing the vorticity field.
- `delT`: Time step size.
- `fn`: Flag indicating which field to use from `varnw`.
- `Nt`: Number of time steps.
- `omL`: Array of shape `(ni, nj, nk, Nt, 3)` to store the Lagrangian vorticity field.
- `varsn`: Array of shape `(ni, nj, nk, Nt, nr)` to store the Lagrangian variables.
- `volsn`: Array of shape `(ni, nj, nk, Nt, 2)` to store the Lagrangian volumes.

## Example
```julia
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

EUL_to_LAG(ni, nj, nk, nr, grdnw, varnw, delT, fn, Nt, omL, varsn, volsn)
"""

function initialize_fields(ni, nj, nr, Nt, grdnw, varnw)
    omL = zeros(ni, nj, Nt, 3)
    varsn = zeros(ni, nj, Nt, nr)
    volsn = ones(ni, nj, Nt, 2)
    
    omL[:, :, 1, 1:3] .= grdnw[:, :, 1:3]
    varsn[:, :, 1, 1:nr] .= varnw[:, :, 1, 1:nr]
    volsn[:, :, 1, 1] .= volsn[:, :, 1, 2]
    
    return omL, varsn, volsn
end

"""
    calculate_surface_position(omL, delT, varsnL)

This function updates the positions of Lagrangian points based on the current velocity field.

# Arguments
- `omL`: A 5-dimensional array representing the Lagrangian point positions.
- `delT`: The time step used for the update.
- `varsnL`: A 4-dimensional array representing the current velocity field.

# Returns
- `omLL`: A 4-dimensional array containing the updated Lagrangian point positions.

# Details
The `calculate_surface_position` function calculates the new positions of the Lagrangian points
after a time step (`delT`). It adds the displacement obtained from the current velocity field (`varsnL`)
to the Lagrangian point positions (`omL`). The updated positions are stored in a new array `omLL`.

"""
function calculate_surface_position(omL, delT, varsnL)
    # Create a copy of the current Lagrangian point positions (omL) and store it in omLL.
    omLL = copy(omL)
    # Multiply the velocity field (varsnL) by the time step (delT) and add the result to
    # the corresponding positions in omLL. This step accounts for the displacement of the
    # Lagrangian points due to the velocity.
    omLL[:, :, 1:3] .+= delT .* varsnL[:, :, 1:3]
    return omLL
end

"""
    calculate_velocity(ni, nj, omLL, grdnw, varnw, fn, iNt)

The calculate_velocity function calculates the velocity values at Lagrangian points based
on their positions and the provided grid information. It also determines the corresponding
values of other variables and volume information associated with the Lagrangian points.

# Arguments
- `ni`: An integer specifying the number of grid points along the x-axis.
- `nj`: An integer specifying the number of grid points along the y-axis.
- `omLL`: A 3-dimensional array representing the Lagrangian point positions (x, y, z).
- `grdnw`: A 4-dimensional array representing the grid information (x, y, z) used for interpolation.
- `varnw`: A 5-dimensional array representing the variable values associated with the grid points.
- `fn`: An integer indicating the selection criteria for variable values.
- `iNt`: An integer specifying the time index.

# Returns
- `varsnL`: A 3-dimensional array representing the variable values at Lagrangian points.
- `volsnL`: A 2-dimensional array representing the volume information at Lagrangian points.

# Details
The calculate_velocity function iterates over each Lagrangian point specified by omLL and
determines the corresponding grid cell indices based on the provided grid information grdnw.
Using these grid cell indices, the function assigns the appropriate variable values from
varnw to varsnL based on the selection criteria defined by fn and the time index iNt.
Similarly, the function assigns the corresponding volume values from volsn to volsnL.
The calculated velocity values, variable values (varsnL), and volume information (volsnL)
are returned as output.

Note: The function assumes a 2D scenario, where only the x and y coordinates are considered
for interpolation and calculation.
"""


function calculate_velocity(ni, nj, omLL, grdnw, varnw, fn, iNt)
    varsnL = zeros(ni, nj, nr)
    volsnL = zeros(ni, nj)
    
    for ii in 1:ni, jj in 1:nj
        for i = 1:ni-1
            if omLL[ii, jj, 1] ≥ grdnw[i, j, 1] && omLL[ii, jj, 1] < grdnw[i+1, j, 1]
                break
            end
        end
        for j = 1:nj-1
            if omLL[ii, jj, 2] ≥ grdnw[i, j, 2] && omLL[ii, jj, 2] < grdnw[i, j+1, 2]
                break
            end
        end
        if fn == 1
            varsnL[ii, jj, 1:nr] = varnw[i, j, 1, 1:nr]
            volsnL[ii, jj] = volsn[i, j, 1, 2]
        else
            varsnL[ii, jj, 1:nr] = varnw[i, j, iNt, 1:nr]
            volsnL[ii, jj] = volsn[i, j, iNt, 2]
        end
    end
    
    return varsnL, volsnL
end

function store_fields(omL, varsn, volsn, omLL, varsnL, volsnL, iNt)
    omL[:, :, iNt, 1:3] .= omLL[:, :, 1:3]
    varsn[:, :, iNt, 1:nr] .= varsnL[:, :, 1:nr]
    volsn[:, :, iNt, 1] .= volsnL[:, :] ./ sum(volsnL[:, :])
end

function EUL_to_LAG(ni, nj, nr, fn, Nt, delT, grdnw, varnw)
    omL, varsn, volsn = initialize_fields(ni, nj, nr, fn, Nt, grdnw, varnw)
    omLL = copy(omL)
    varsnL = copy(varsn)
    volsnL = copy(volsn)
    
    for iNt = 2:Nt
        omLL = calculate_surface_position(omLL, delT, varsnL)
        varsnL, volsnL = calculate_velocity(ni, nj, omLL, grdnw, varnw, fn, iNt)
        store_fields(omL, varsn, volsn, omLL, varsnL, volsnL, iNt)
        
        println(iNt, "/", Nt)
    end
    
    return omL, varsn, volsn
end

ni = 10
nj = 10
nr = 5
Nt = 3
fn = 2
delT = 0.1

grdnw = rand(ni, nj, 3)
varnw = rand(ni, nj, fn, nr)
omLL, varsnL, volsnL = initialize_fields(ni, nj, nr, Nt, grdnw, varnw)
omLL = calculate_surface_position(omLL, delT, varsnL)
varsnL, volsnL = calculate_velocity(ni, nj, omLL, grdnw, varnw, fn, iNt)


# omL, varsn, volsn = EUL_to_LAG(ni, nj, nr, fn, Nt, delT, grdnw, varnw)
