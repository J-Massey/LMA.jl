using LinearAlgebra

function lagrangian_POD_analysis(trajectory_matrix)
    dim, nx, ny, nt = size(trajectory_matrix)
    δt = 1.0  # Time step

    # Convert Eulerian snapshots to Lagrangian snapshots
    δτ = δt  # Forward passage of time
    χ = zeros(dim, nx, ny, nt)
    τ = zeros(nt)
    U = copy(trajectory_matrix)  # Eulerian snapshots

    for i in 1:nt
        U[:, :, :, i] = trajectory_matrix[:, :, :, i]
        χ[:, :, :, i] = δτ * U[:, :, :, i]
        τ[i] = i * δτ
    end

    # Lagrangian POD
    C = zeros(nt, nt)
    W = ones(nt, nt)  # Weight matrix, assuming unity weights

    for i in 1:nt
        for j in 1:nt
            C[i, j] = dot(vec(U[:, :, :, i]), W[i, j] * vec(U[:, :, :, j]))
        end
    end

    λ, Ψ = eigen(C)
    Φ = zeros(dim, nx, ny, nt)

    for n in 1:nt
        Φ[:, :, :, n] = sqrt(λ[n]) * sum(U[:, :, :, i] .* Ψ[i, n] for i in 1:nt)
    end

    return Φ
end

# Usage example
nx = 10  # Number of grid points in x-direction
ny = 10  # Number of grid points in y-direction
nt = 100  # Number of snapshots
dim = 2  # Number of dimensions (x and y)
# trajectory_matrix = rand(dim, nx, ny, nt)  # Sample trajectory matrix
traj_x = npzread("data/data/lagrangian_trajectories_x.npy")
traj_y = npzread("data/data/lagrangian_trajectories_y.npy")

# Combine traj_x and traj_x to get trajectory_matrix
trajectory_matrix = cat(traj_x, traj_y, dims = 4)
trajectory_matrix = permutedims(trajectory_matrix, [4, 1, 2, 3])

modes = lagrangian_POD_analysis(trajectory_matrix)
npzwrite("data/data/lagrangian_POD_modes.npy", modes)