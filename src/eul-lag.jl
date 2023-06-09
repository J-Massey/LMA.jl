using Interpolations
using BenchmarkTools
using GLMakie
using DynamicalSystems

"""
    eulerian_to_lagrangian(flow_snapshots, particle_positions)

Converts flow snapshots from an Eulerian frame to a Lagrangian frame.

# Arguments
- `flow_snapshots`: Array containing flow snapshots at different time steps.
- `particle_positions`: Array containing initial positions of fluid particles.

# Returns
- `lagrangian_trajectories`: Array containing Lagrangian trajectories of fluid particles.

"""
# Test data
xlims, ylims = (-0.35, 2), (-0.35, 0.35)
nx,ny,nt = 1080,1920,10

pxs = LinRange(xlims..., nx)
pys = LinRange(ylims..., ny)

flow_snapshots_u = rand(nx, ny, nt)
flow_snapshots_v = rand(nx, ny, nt)

# Now set up a uniform grid of particles
xgrid = repeat(pxs, 1, size(flow_snapshots_u, 2))
ygrid = repeat(pys, 1, size(flow_snapshots_u, 1))'

# For the initial time step, the trajectory is the same as the position
lagrangian_trajectories_x = repeat(xgrid, 1, 1, nt)
lagrangian_trajectories_y = repeat(ygrid, 1, 1, nt)

# Iterate over time steps
for t in 2:nt
    # Extract velocity field for current snapshot
    U =  @view flow_snapshots_u[:,:,t]
    space_interpolation_u = interpolate((pxs, pys), U, Gridded(Linear()))
    V = @view flow_snapshots_v[:,:,t]
    space_interpolation_y = interpolate((pxs, pys), V, Gridded(Linear()))

    for p1 in 1:nx
        for p2 in 1:ny
            x = pxs[p1]
            y = pys[p2]
            ui = space_interpolation_u(x, y)
            vi = space_interpolation_y(x, y)
            prev_x = lagrangian_trajectories_x[p1, p2, t - 1]
            prev_y = lagrangian_trajectories_y[p1, p2, t - 1]
            lagrangian_trajectories_x[p1, p2, t] = prev_x + ui # ui*dx, so positional
            lagrangian_trajectories_y[p1, p2, t] = prev_y + vi
        end
    end
end

for ti in 1:nt
    fig = Figure()
    ax = Axis(fig[1, 1])
    xlims!(ax, xlims...)
    ylims!(ax, ylims...)
    println("Plotting time step $ti")
    scatter!(ax, vec(lagrangian_trajectories_x[:,:,ti]), vec(lagrangian_trajectories_y[:,:,ti]), markersize = 0.5, color = :blue)
    save("figures/lagrangian_trajectories_$ti.png", fig)
end
