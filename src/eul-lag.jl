using Interpolations
using GLMakie

"""
    eulerian_to_lagrangian(flow_snapshots, particle_positions)

Converts flow snapshots from an Eulerian frame to a Lagrangian frame.

# Arguments
- `flow_snapshots`: Array containing flow snapshots at different time steps.
- `particle_positions`: Array containing initial positions of fluid particles.

# Returns
- `lagrangian_trajectories`: Array containing Lagrangian trajectories of fluid particles.

"""
# Example data, start 1D
n,nt = 100,10
flow_snapshots_u = rand(n, nt)
# Now set up a uniform grid of particles
px = [Float64(i) for i in 1:n]

# function eulerian_to_lagrangian(flow_snapshots, particle_positions)
num_snapshots = nt
num_particles = n

# For the initial time step, the trajectory is the same as the position
lagrangian_trajectories_x = repeat(px, 1, num_snapshots)

# Iterate over time steps
for t in 2:num_snapshots
    # Extract velocity field for current snapshot
    velocity_field = @view flow_snapshots_u[:,t]
    space_interpolation = interpolate(velocity_field, BSpline(Linear()))
    for p1 in 1:n
        x = px[p1]
        u = space_interpolation(x)
        prev_x = lagrangian_trajectories_x[p1, t - 1]
        lagrangian_trajectories_x[p1, t] = x + u
    end
end

println(size(lagrangian_trajectories_x))
