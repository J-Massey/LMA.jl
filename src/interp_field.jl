using Interpolations

# Example data
velocity_field = rand(100, 100)  # Velocity field matrix of size 100x100x3
positions = rand(100, 100)  # Array of positions of size 100x100

# Define interpolation objects for each component of the velocity field
interpolations = interpolate(velocity_field[:, :], BSpline(Linear()))



# Example usage
x_pos = 50  # X-position to interpolate
y_pos = 75  # Y-position to interpolate
interpolated_velocities = interpolations(x_pos, y_pos)

interpolated_velocities = interpolate_velocity((x_pos, y_pos))

# Access the interpolated velocities
vx, vy, vz = interpolated_velocities
println("Interpolated velocities at position ($x_pos, $y_pos):")
println("Vx: $vx")
println("Vy: $vy")
println("Vz: $vz")
