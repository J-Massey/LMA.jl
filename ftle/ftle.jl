using DynamicalSystems
using CairoMakie


function flow_field(flow_field, n) # here `n` is "time", but we don't use it.
    xn, yn = flow_field[n, :, :] # system state
    return SVector(xn, yn)
end

# Example data
n,m,t = 100,100,10
flow = rand(n, m, t)
particle_positions = [Float64[i, j] for i in 1:n, j in 1:m]

henon = DeterministicIteratedMap(flow_field, flow[1, :, :], particle_positions)

total_time = 10_000
X, t = trajectory(henon, total_time)

scatter(X[:, 1], X[:, 2])


# ds_struct = convert_flow_field_to_trajectory(flow_field_data)
