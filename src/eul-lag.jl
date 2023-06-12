using Interpolations
# using BenchmarkTools
using GLMakie
# using DynamicalSystems
using NPZ
using Statistics

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

# Load flow snapshots
flow_snapshots_u = npzread("data/data/u.npy")
flow_snapshots_u = permutedims(flow_snapshots_u, [3, 2, 1])
flow_snapshots_v = npzread("data/data/v.npy")
flow_snapshots_v = permutedims(flow_snapshots_v, [3, 2, 1])
# take away the time mean
flow_snapshots_u = flow_snapshots_u .- mean(flow_snapshots_u, dims = 3)
flow_snapshots_v = flow_snapshots_v .- mean(flow_snapshots_v, dims = 3)

function testplot()
    # Quick test plot
    f = Figure()
    ax = Axis(f[1, 1])
    co = contourf!(ax,
        collect(flow_snapshots_u[:, :, 1]),
        xlabel=L"x", ylabel=L"y", title=L"ω_z",
        # levels = range(-0.1, 1.5, length = 44),
        # colormap=:icefire,
        extendlow = :auto, extendhigh = :auto,
    )
    tightlimits!(ax)
    return f
end

nx, ny, nt = size(flow_snapshots_u)
dt = 4/nt

pxs = LinRange(xlims..., nx)
pys = LinRange(ylims..., ny)

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
            lagrangian_trajectories_x[p1, p2, t] = prev_x + ui*dt
            lagrangian_trajectories_y[p1, p2, t] = prev_y + vi*dt
        end
    end
end

for ti in 1:nt
    fig = Figure()
    ax = Axis(fig[1, 1])
    xlims!(ax, -0.35, 2.5)
    ylims!(ax, -0.35, 0.35)
    println("Plotting time step $ti")
    scatter!(ax, vec(lagrangian_trajectories_x[:,:,ti]), vec(lagrangian_trajectories_y[:,:,ti]), markersize = 0.5, color = :blue)
    save("figures/lagrangian_trajectories_$ti.png", fig)
end
