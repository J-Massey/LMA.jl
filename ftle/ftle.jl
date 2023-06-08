using DynamicalSystems
using CairoMakie


function flow_field(flow_field, n) # here `n` is "time", but we don't use it.
    xn, yn = flow_field[n, :, :] # system state
    return SVector(xn, yn)
end

u0 = [0.2, 0.3]
p0 = [1.4, 0.3]

t = 20; nx=1000; ny=1000
flow_field_size = (t, nx, ny)
flow_field_data = rand(t, nx, ny) 

 # Your flow field data goes here

henon = DeterministicIteratedMap(flow_field, flow_field[1, :, :])

total_time = 10_000
X, t = trajectory(henon, total_time)

scatter(X[:, 1], X[:, 2])


# ds_struct = convert_flow_field_to_trajectory(flow_field_data)
