"""
This file includes expoential integrators for states and unitaries
"""

using ExponentialAction

function exp_eigen(G::AbstractMatrix)
    Ĥ = Hermitian(Matrix(Isomorphisms.H(G)))
    λ, V = eigen(Ĥ)
    expG = Isomorphisms.iso(sparse(V * Diagonal(exp.(-im * λ)) * V'))
    droptol!(expG, 1e-12)
    return expG
end

# ----------------------------------------------------------------------------- #
#                         Quantum Exponential Integrators                       #
# ----------------------------------------------------------------------------- #

abstract type QuantumExponentialIntegrator <: QuantumIntegrator end

# ----------------------------------------------------------------------------- #
#                         Unitary Exponential Integrator                        #
# ----------------------------------------------------------------------------- #

struct UnitaryExponentialIntegrator <: QuantumExponentialIntegrator
    unitary_components::Vector{Int}
    drive_components::Vector{Int}
    timestep::Union{Real, Int}
    freetime::Bool
    n_drives::Int
    ketdim::Int
    dim::Int
    zdim::Int
    autodiff::Bool
    G::Function

    function UnitaryExponentialIntegrator(
        unitary_name::Symbol,
        drive_name::Union{Symbol, Tuple{Vararg{Symbol}}},
        sys::QuantumSystem,
        traj::NamedTrajectory;
        autodiff::Bool=false
    )
        dim = traj.dims[unitary_name]

        ketdim = Int(sqrt(dim ÷ 2))

        unitary_components = traj.components[unitary_name]

        if drive_name isa Tuple
            drive_components = vcat((traj.components[s] for s ∈ drive_name)...)
        else
            drive_components = traj.components[drive_name]
        end

        n_drives = length(drive_components)

        @assert all(diff(drive_components) .== 1) "controls must be in order"

        freetime = traj.timestep isa Symbol

        if freetime
            timestep = traj.components[traj.timestep][1]
        else
            timestep = traj.timestep
        end

        return new(
            unitary_components,
            drive_components,
            timestep,
            freetime,
            n_drives,
            ketdim,
            dim,
            traj.dim,
            autodiff,
            sys.G
        )
    end
end

function (integrator::UnitaryExponentialIntegrator)(
    traj::NamedTrajectory;
    unitary_name::Union{Nothing, Symbol}=nothing,
    drive_name::Union{Nothing, Symbol, Tuple{Vararg{Symbol}}}=nothing,
    G::Function=integrator.G,
    autodiff::Bool=integrator.autodiff
)
    @assert !isnothing(unitary_name) "unitary_name must be provided"
    @assert !isnothing(drive_name) "drive_name must be provided"
    return UnitaryExponentialIntegrator(
        unitary_name,
        drive_name,
        G,
        traj;
        autodiff=autodiff
    )
end

function get_comps(P::UnitaryExponentialIntegrator, traj::NamedTrajectory)
    if P.freetime
        return P.unitary_components, P.drive_components, traj.components[traj.timestep]
    else
        return P.unitary_components, P.drive_components
    end
end

# ------------------------------ Integrator --------------------------------- #

@views function (ℰ::UnitaryExponentialIntegrator)(
    zₜ::AbstractVector,
    zₜ₊₁::AbstractVector,
    t::Int
)
    Ũ⃗ₜ₊₁ = zₜ₊₁[ℰ.unitary_components]
    Ũ⃗ₜ = zₜ[ℰ.unitary_components]
    aₜ = zₜ[ℰ.drive_components]

    if ℰ.freetime
        Δtₜ = zₜ[ℰ.timestep]
    else
        Δtₜ = ℰ.timestep
    end

    # return Ũ⃗ₜ₊₁ - expv(Δtₜ, I(ℰ.ketdim) ⊗ ℰ.G(aₜ), Ũ⃗ₜ)
    return Ũ⃗ₜ₊₁ - (I(ℰ.ketdim) ⊗ exp_eigen(Δtₜ * ℰ.G(aₜ))) * Ũ⃗ₜ
end

@views function jacobian(
    ℰ::UnitaryExponentialIntegrator,
    zₜ::AbstractVector,
    zₜ₊₁::AbstractVector,
    t::Int
)
    # get the state and control vectors
    Ũ⃗ₜ = zₜ[ℰ.unitary_components]
    aₜ = zₜ[ℰ.drive_components]

    # obtain the timestep
    if ℰ.freetime
        Δtₜ = zₜ[ℰ.timestep]
    else
        Δtₜ = ℰ.timestep
    end

    # compute the generator
    Gₜ = ℰ.G(aₜ)

    Id = I(ℰ.ketdim)

    expGₜ = exp_eigen(Δtₜ * Gₜ)

    ∂Ũ⃗ₜ₊₁ℰ = sparse(I, ℰ.dim, ℰ.dim)
    ∂Ũ⃗ₜℰ = -Id ⊗ expGₜ

    ∂aₜℰ = ForwardDiff.jacobian(
        a -> -expv(Δtₜ, Id ⊗ ℰ.G(a), Ũ⃗ₜ),
        aₜ
    )

    if ℰ.freetime
        ∂Δtₜℰ = -(Id ⊗ (Gₜ * expGₜ)) * Ũ⃗ₜ
        return ∂Ũ⃗ₜℰ, ∂Ũ⃗ₜ₊₁ℰ, ∂aₜℰ, ∂Δtₜℰ
    else
        return ∂Ũ⃗ₜℰ, ∂Ũ⃗ₜ₊₁ℰ, ∂aₜℰ
    end
end

struct QuantumStateExponentialIntegrator <: QuantumExponentialIntegrator
    state_components::Vector{Int}
    drive_components::Vector{Int}
    timestep::Union{Real, Int}
    freetime::Bool
    n_drives::Int
    ketdim::Int
    dim::Int
    zdim::Int
    autodiff::Bool
    G::Function

    function QuantumStateExponentialIntegrator(
        state_name::Symbol,
        drive_name::Union{Symbol, Tuple{Vararg{Symbol}}},
        sys::QuantumSystem,
        traj::NamedTrajectory;
        autodiff::Bool=false
    )
        dim = traj.dims[state_name]
        ketdim = dim ÷ 2

        state_components = traj.components[state_name]

        if drive_name isa Tuple
            drive_components = vcat((traj.components[s] for s ∈ drive_name)...)
        else
            drive_components = traj.components[drive_name]
        end

        n_drives = length(drive_components)

        @assert all(diff(drive_components) .== 1) "controls must be in order"

        freetime = traj.timestep isa Symbol

        if freetime
            timestep = traj.components[traj.timestep][1]
        else
            timestep = traj.timestep
        end

        return new(
            state_components,
            drive_components,
            timestep,
            freetime,
            n_drives,
            ketdim,
            dim,
            traj.dim,
            autodiff,
            sys.G
        )
    end
end

function get_comps(P::QuantumStateExponentialIntegrator, traj::NamedTrajectory)
    if P.freetime
        return P.state_components, P.drive_components, traj.components[traj.timestep]
    else
        return P.state_components, P.drive_components
    end
end

# function (integrator::QuantumStateExponentialIntegrator)(
#     traj::NamedTrajectory;
#     state_name::Union{Nothing, Symbol}=nothing,
#     drive_name::Union{Nothing, Symbol, Tuple{Vararg{Symbol}}}=nothing,
#     G::Function=integrator.G,
#     autodiff::Bool=integrator.autodiff
# )
#     @assert !isnothing(state_name) "state_name must be provided"
#     @assert !isnothing(drive_name) "drive_name must be provided"
#     return QuantumStateExponentialIntegrator(
#         state_name,
#         drive_name,
#         G,
#         traj;
#         autodiff=autodiff
#     )
# end

# ------------------------------ Integrator --------------------------------- #

@views function (ℰ::QuantumStateExponentialIntegrator)(
    zₜ::AbstractVector,
    zₜ₊₁::AbstractVector,
    t::Int
)
    ψ̃ₜ₊₁ = zₜ₊₁[ℰ.state_components]
    ψ̃ₜ = zₜ[ℰ.state_components]
    aₜ = zₜ[ℰ.drive_components]

    if ℰ.freetime
        Δtₜ = zₜ[ℰ.timestep]
    else
        Δtₜ = ℰ.timestep
    end

    return ψ̃ₜ₊₁ - expv(Δtₜ, ℰ.G(aₜ), ψ̃ₜ)
end

@views function jacobian(
    ℰ::QuantumStateExponentialIntegrator,
    zₜ::AbstractVector,
    zₜ₊₁::AbstractVector,
    t::Int
)
    # get the state and control vectors
    ψ̃ₜ = zₜ[ℰ.state_components]
    aₜ = zₜ[ℰ.drive_components]

    # obtain the timestep
    if ℰ.freetime
        Δtₜ = zₜ[ℰ.timestep]
    else
        Δtₜ = ℰ.timestep
    end

    # compute the generator
    Gₜ = ℰ.G(aₜ)

    expGₜ = sparse(exp_eigen(Δtₜ * Gₜ))

    ∂ψ̃ₜ₊₁ℰ = sparse(I, ℰ.dim, ℰ.dim)
    ∂ψ̃ₜℰ = -expGₜ

    ∂aₜℰ = ForwardDiff.jacobian(
        a -> -expv(Δtₜ, ℰ.G(a), ψ̃ₜ),
        aₜ
    )

    if ℰ.freetime
        ∂Δtₜℰ = -Gₜ * (expGₜ * ψ̃ₜ)
        return ∂ψ̃ₜℰ, ∂ψ̃ₜ₊₁ℰ, ∂aₜℰ, ∂Δtₜℰ
    else
        return ∂ψ̃ₜℰ, ∂ψ̃ₜ₊₁ℰ, ∂aₜℰ
    end
end

# ----------------------------------------------------------------------------- #
#                Density Operator Exponential Integrator                        #
# ----------------------------------------------------------------------------- #

struct DensityOperatorExponentialIntegrator <: QuantumExponentialIntegrator
    density_operator_components::Vector{Int}
    drive_components::Vector{Int}
    timestep::Union{Real, Int}
    freetime::Bool
    n_drives::Int
    ketdim::Int
    dim::Int
    zdim::Int
    autodiff::Bool
    G::Function

    function DensityOperatorExponentialIntegrator(
        density_operator_name::Symbol,
        drive_name::Union{Symbol, Tuple{Vararg{Symbol}}},
        sys::AbstractQuantumSystem,
        traj::NamedTrajectory;
        autodiff::Bool=false
    )
        dim = traj.dims[density_operator_name]
        ketdim = size(sys.H(zeros(sys.n_drives)), 1)

        density_operator_components = traj.components[density_operator_name]

        println(length(density_operator_components))

        if drive_name isa Tuple
            drive_components = vcat((traj.components[s] for s ∈ drive_name)...)
        else
            drive_components = traj.components[drive_name]
        end

        n_drives = length(drive_components)

        @assert all(diff(drive_components) .== 1) "controls must be in order"

        freetime = traj.timestep isa Symbol

        if freetime
            timestep = traj.components[traj.timestep][1]
        else
            timestep = traj.timestep
        end

        return new(
            density_operator_components,
            drive_components,
            timestep,
            freetime,
            n_drives,
            ketdim,
            dim,
            traj.dim,
            autodiff,
            sys.G
        )
    end
end

# ------------------------------ Integrator --------------------------------- #

@views function (ℒ::DensityOperatorExponentialIntegrator)(
    zₜ::AbstractVector,
    zₜ₊₁::AbstractVector,
    t::Int
)
    ρ⃗̃ₜ₊₁ = zₜ₊₁[ℒ.density_operator_components]
    ρ⃗̃ₜ = zₜ[ℒ.density_operator_components]
    aₜ = zₜ[ℒ.drive_components]

    if ℒ.freetime
        Δtₜ = zₜ[ℒ.timestep]
    else
        Δtₜ = ℒ.timestep
    end

    return ρ⃗̃ₜ₊₁ - expv(Δtₜ, ℒ.G(aₜ), ρ⃗̃ₜ)
end

# ------------------------------ Jacobian --------------------------------- #

@views function jacobian(
    ℒ::DensityOperatorExponentialIntegrator,
    zₜ::AbstractVector,
    zₜ₊₁::AbstractVector,
    t::Int
)
    # get the state and control vectors
    ρ⃗̃ₜ = zₜ[ℒ.density_operator_components]
    aₜ = zₜ[ℒ.drive_components]

    # obtain the timestep
    if ℒ.freetime
        Δtₜ = zₜ[ℒ.timestep]
    else
        Δtₜ = ℒ.timestep
    end

    # compute the generator
    Gₜ = ℒ.G(aₜ)

    # compute the exponential
    expGₜ = sparse(exp(Δtₜ * Matrix(Gₜ)))

    ∂ρ⃗̃ₜ₊₁ℒ = sparse(I, ℒ.dim, ℒ.dim)

    ∂ρ⃗̃ₜℒ = -expGₜ

    ∂aₜℒ = ForwardDiff.jacobian(
        a -> -expv(Δtₜ, ℒ.G(a), ρ⃗̃ₜ),
        aₜ
    )

    if ℒ.freetime
        ∂Δtₜℒ = -Gₜ * (expGₜ * ρ⃗̃ₜ)
        return ∂ρ⃗̃ₜℒ, ∂ρ⃗̃ₜ₊₁ℒ, ∂aₜℒ, ∂Δtₜℒ
    else
        return ∂ρ⃗̃ₜℒ, ∂ρ⃗̃ₜ₊₁ℒ, ∂aₜℒ
    end
end

function get_comps(P::DensityOperatorExponentialIntegrator, traj::NamedTrajectory)
    if P.freetime
        return P.density_operator_components, P.drive_components, traj.components[traj.timestep]
    else
        return P.density_operator_components, P.drive_components
    end
end

# ******************************************************************************* #

@testitem "testing UnitaryExponentialIntegrator" begin
    using NamedTrajectories
    using PiccoloQuantumObjects
    using FiniteDiff

    T = 100
    H_drift = GATES[:Z]
    H_drives = [GATES[:X], GATES[:Y]]
    n_drives = length(H_drives)

    sys = QuantumSystem(H_drift, H_drives)

    U_init = GATES[:I]
    U_goal = GATES[:X]

    Ũ⃗_init = operator_to_iso_vec(U_init)
    Ũ⃗_goal = operator_to_iso_vec(U_goal)

    dt = 0.1


    Z = NamedTrajectory(
        (
            # Ũ⃗ = unitary_geodesic(U_goal, T),
            Ũ⃗ = randn(length(Ũ⃗_init), T),
            a = randn(n_drives, T),
            da = randn(n_drives, T),
            Δt = fill(dt, 1, T),
        ),
        controls=(:da,),
        timestep=:Δt,
        goal=(Ũ⃗ = Ũ⃗_goal,)
    )

    ℰ = UnitaryExponentialIntegrator(:Ũ⃗, :a, sys, Z)


    ∂Ũ⃗ₜℰ, ∂Ũ⃗ₜ₊₁ℰ, ∂aₜℰ, ∂Δtₜℰ = jacobian(ℰ, Z[1].data, Z[2].data, 1)

    ∂ℰ_finitediff= FiniteDiff.finite_difference_jacobian(
        zz -> ℰ(zz[1:Z.dim], zz[Z.dim+1:end], 1),
        [Z[1].data; Z[2].data]
    )

    @test isapprox(∂Ũ⃗ₜℰ, ∂ℰ_finitediff[:,1:ℰ.dim]; atol=1e-6)
    @test isapprox(∂Ũ⃗ₜ₊₁ℰ, ∂ℰ_finitediff[:,Z.dim .+ (1:ℰ.dim)]; atol=1e-6)
    @test isapprox(∂aₜℰ, ∂ℰ_finitediff[:,Z.components.a]; atol=1e-6)
    @test isapprox(∂Δtₜℰ, ∂ℰ_finitediff[:,Z.components.Δt]; atol=1e-6)
end

@testitem "testing QuantumStateExponentialIntegrator" begin
    using NamedTrajectories
    using PiccoloQuantumObjects
    using ForwardDiff

    T = 100
    H_drift = GATES[:Z]
    H_drives = [GATES[:X], GATES[:Y]]
    n_drives = length(H_drives)

    sys = QuantumSystem(H_drift, H_drives)

    U_init = GATES[:I]
    U_goal = GATES[:X]

    ψ̃_init = ket_to_iso([1.0, 0.0])
    ψ̃_goal = ket_to_iso([0.0, 1.0])

    dt = 0.1

    Z = NamedTrajectory(
        (
            # ψ̃ = linear_interpolation(ψ̃_init, ψ̃_goal, T),
            ψ̃ = randn(length(ψ̃_init), T),
            a = randn(n_drives, T),
            da = randn(n_drives, T),
            Δt = fill(dt, 1, T),
        ),
        controls=(:da,),
        timestep=:Δt,
        goal=(ψ̃ = ψ̃_goal,)
    )

    ℰ = QuantumStateExponentialIntegrator(:ψ̃, :a, sys, Z)

    ∂ψ̃ₜℰ, ∂ψ̃ₜ₊₁ℰ, ∂aₜℰ, ∂Δtₜℰ = jacobian(ℰ, Z[1].data, Z[2].data, 1)

    ∂ℰ_forwarddiff = ForwardDiff.jacobian(
        zz -> ℰ(zz[1:Z.dim], zz[Z.dim+1:end], 1),
        [Z[1].data; Z[2].data]
    )

    @test ∂ψ̃ₜℰ ≈ ∂ℰ_forwarddiff[:, 1:ℰ.dim]
    @test ∂ψ̃ₜ₊₁ℰ ≈ ∂ℰ_forwarddiff[:, Z.dim .+ (1:ℰ.dim)]
    @test ∂aₜℰ ≈ ∂ℰ_forwarddiff[:, Z.components.a]
    @test ∂Δtₜℰ ≈ ∂ℰ_forwarddiff[:, Z.components.Δt]
end

@testitem "testing DensityOperatorExponentialIntegrator" begin
    using NamedTrajectories
    using PiccoloQuantumObjects
    using ForwardDiff

    T = 100
    H_drift = GATES[:Z]
    H_drives = [GATES[:X], GATES[:Y]]
    n_drives = length(H_drives)

    ψ0 = [1.0, 0.0]
    ψ1 = [0.0, 1.0]

    sys = OpenQuantumSystem(H_drift, H_drives, [ψ0 * ψ1'])


    ρ_init = ψ0 * ψ0'
    ρ_goal = ψ1 * ψ1'

    ρ⃗̃_init = density_to_iso_vec(ρ_init)
    ρ⃗̃_goal = density_to_iso_vec(ρ_goal)

    dt = 0.1

    Z = NamedTrajectory(
        (
            # ρ⃗̃ = linear_interpolation(ρ⃗̃_init, ρ⃗̃_goal, T),
            ρ⃗̃ = randn(length(ρ⃗̃_init), T),
            a = randn(n_drives, T),
            da = randn(n_drives, T),
            Δt = fill(dt, 1, T),
        ),
        controls=(:da, :Δt),
        timestep=:Δt,
        goal=(ρ⃗̃ = ρ⃗̃_goal,)
    )

    ℒ = DensityOperatorExponentialIntegrator(:ρ⃗̃, :a, sys, Z)

    ∂ρ⃗̃ₜℒ, ∂ρ⃗̃ₜ₊₁ℒ, ∂aₜℒ, ∂Δtₜℒ = jacobian(ℒ, Z[1].data, Z[2].data, 1)

    ∂ℒ_forwarddiff = ForwardDiff.jacobian(
        zz -> ℒ(zz[1:Z.dim], zz[Z.dim+1:end], 1),
        [Z[1].data; Z[2].data]
    )

    @test ∂ρ⃗̃ₜℒ ≈ ∂ℒ_forwarddiff[:, 1:ℒ.dim]
    @test ∂ρ⃗̃ₜ₊₁ℒ ≈ ∂ℒ_forwarddiff[:, Z.dim .+ (1:ℒ.dim)]
    @test ∂aₜℒ ≈ ∂ℒ_forwarddiff[:, Z.components.a]
    @test ∂Δtₜℒ ≈ ∂ℒ_forwarddiff[:, Z.components.Δt]
end
