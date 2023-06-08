using Interpolations

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
n,m,t = 100,100,10
flow_snapshots_u = rand(n, m, t)
# Now set up a uniform grid of particles
px = [i for i in 1:n]
py = [j for j in 1:m]

# function eulerian_to_lagrangian(flow_snapshots, particle_positions)
num_snapshots = t
num_particles = n*m
lagrangian_trajectories = copy(flow_snapshots_u)

# Iterate over time steps
for t in 1:num_snapshots
    # Extract velocity field for current snapshot
    velocity_field = @view flow_snapshots_u[:,:,t]
    space_interpolation = interpolate(velocity_field, BSpline(Linear()))
    for p1 in 1:n
        for p2 in 1:m
            x = px[p1]
            y = py[p2]
            u = space_interpolation(x, y)
            # println(lagrangian_trajectories[p1, p2, t]) = x + u
        end
    end
    
    # Iterate over particles
    # for p in 1:num_particles
    #     # Get particle position in the Eulerian frame
    #     x, y = particle_positions[p, :]
        
    #     # Integrate velocity to obtain particle trajectory
    #     if t == 1
    #         # For the initial time step, the trajectory is the same as the position
    #         lagrangian_trajectories[p, t] = (x, y)
    #     else
    #         # Integrate using the Euler method
    #         dt = t  # Assuming time step is 1 unit
    #         prev_x, prev_y, prev_z = lagrangian_trajectories[p, t - 1]
    #         lagrangian_trajectories[p, t] = (prev_x + dt * u, prev_y + dt * vy, prev_z + dt * vz)
    #     end
    # end
end
    
#     return lagrangian_trajectories
# end

"""
    main()

Example usage of `eulerian_to_lagrangian` function.
"""
# function main()
# Example data
flow_snapshots = rand(100, 100, 10)
# Now set up a uniform grid of particles
particle_positions = [
    (i, j, k) for i in 1:size(flow_snapshots, 1),
    j in 1:size(flow_snapshots, 2),
    k in 1:size(flow_snapshots, 3)
    ]

# Call the function to convert to Lagrangian frame
lagrangian_trajectories = eulerian_to_lagrangian(flow_snapshots, particle_positions)

# Access the particle trajectories at specific time step
t = 3
particle_1_trajectory = lagrangian_trajectories[1, t]
println("Particle 1 trajectory at time step $t: $particle_1_trajectory")
# end
