module Dynamics

export AbstractDynamics
export QuantumDynamics

export dynamics
export dynamics_jacobian
export dynamics_components

using ..Integrators

using TrajectoryIndexingUtils
using NamedTrajectories
using LinearAlgebra
using SparseArrays
using ForwardDiff


abstract type AbstractDynamics end

"""
    QuantumDynamics <: AbstractDynamics
"""
struct QuantumDynamics <: AbstractDynamics
    integrators::Union{Nothing, Vector{<:AbstractIntegrator}}
    F::Function
    ∂F::Function
    ∂F_structure::Vector{Tuple{Int, Int}}
    μ∂²F::Union{Function, Nothing}
    μ∂²F_structure::Vector{Tuple{Int, Int}}
    dim::Int
end

function NullQuantumDynamics()
    return QuantumDynamics(
        nothing,
        _ -> nothing,
        _ -> nothing,
        [],
        nothing,
        [],
        0
    )
end

function dynamics_components(integrators::Vector{<:AbstractIntegrator})
    dynamics_comps = []
    comp_mark = 0
    for integrator ∈ integrators
        integrator_comps = (comp_mark + 1):(comp_mark + integrator.dim)
        push!(dynamics_comps, integrator_comps)
        comp_mark += integrator.dim
    end
    return dynamics_comps
end

function dynamics(
    integrators::Vector{<:AbstractIntegrator}
)
    dynamics_comps = dynamics_components(integrators)
    dynamics_dim = sum(integrator.dim for integrator ∈ integrators)
    function f(zₜ, zₜ₊₁, t)
        δ = Vector{eltype(zₜ)}(undef, dynamics_dim)
        for (integrator, integrator_comps) ∈ zip(integrators, dynamics_comps)
            δ[integrator_comps] = integrator(zₜ, zₜ₊₁, t)
        end
        return δ
    end
    return f
end

function dynamics_jacobian(
    integrators::Vector{<:AbstractIntegrator}
)
    dynamics_comps = dynamics_components(integrators)
    dynamics_dim = sum(integrator.dim for integrator ∈ integrators)
    zdim = integrators[1].zdim
    @views function ∂f(zₜ, zₜ₊₁, t)
        ∂ = spzeros(eltype(zₜ), dynamics_dim, 2zdim)
        for (integrator, integrator_comps) ∈ zip(integrators, dynamics_comps)
            ∂[integrator_comps, 1:2zdim] =
                Integrators.jacobian(integrator, zₜ, zₜ₊₁, t)
        end
        return ∂
    end
    return ∂f
end


function QuantumDynamics(
    integrators::Vector{<:AbstractIntegrator},
    traj::NamedTrajectory;
    eval_hessian=true,
    verbose=false
)
    if length(integrators) == 0
        if verbose
            println("        constructing Null dynamics function...")
        end

        return NullQuantumDynamics()
    end

    if verbose
        println("        constructing knot point dynamics functions...")
    end

    f = dynamics(integrators)

    ∂f = dynamics_jacobian(integrators)

    if eval_hessian
        error("Hessians not implemented")
    else
        μ∂²f = nothing
    end

    dynamics_dim = sum(integrator.dim for integrator ∈ integrators)

    if eval_hessian
        error("Hessians not implemented")
    else
        ∂f_structure = [(i, j) for i = 1:dynamics_dim, j = 1:2traj.dim]
        ∂F_structure = Vector{Tuple{Int,Int}}(undef, length(∂f_structure) * (traj.T - 1))
        for t = 1:traj.T-1
            ∂fₜ_structure = [
                (
                    i + index(t, 0, dynamics_dim),
                    j + index(t, 0, traj.dim)
                ) for (i, j) ∈ ∂f_structure
            ]
            ∂F_structure[slice(t, length(∂f_structure))] = ∂fₜ_structure
        end

        μ∂²F_structure = []
    end

    if verbose
        println("        constructing full dynamics derivative functions...")
    end

    @views function F(Z⃗::AbstractVector{R}) where R <: Real
        δ = Vector{R}(undef, dynamics_dim * (traj.T - 1))
        Threads.@threads for t = 1:traj.T-1
            zₜ = Z⃗[slice(t, traj.dim)]
            zₜ₊₁ = Z⃗[slice(t + 1, traj.dim)]
            δ[slice(t, dynamics_dim)] = f(zₜ, zₜ₊₁, t)
        end
        return δ
    end

    @views function ∂F(Z⃗::AbstractVector{R}) where R <: Real
        ∂f_nnz = dynamics_dim * 2traj.dim
        ∂s = zeros(R, ∂f_nnz * (traj.T - 1))
        Threads.@threads for t = 1:traj.T-1
            zₜ = Z⃗[slice(t, traj.dim)]
            zₜ₊₁ = Z⃗[slice(t + 1, traj.dim)]
            ∂fₜ = ∂f(zₜ, zₜ₊₁, t)
            for (k, (i, j)) ∈ enumerate(Iterators.product(1:dynamics_dim, 1:2traj.dim))
                ∂s[index(t, k, ∂f_nnz)] = ∂fₜ[i, j]
            end
        end
        return ∂s
    end

    if eval_hessian
        @views μ∂²F = (Z⃗::AbstractVector{<:Real}, μ⃗::AbstractVector{<:Real}) -> begin
            μ∂²s = Vector{eltype(Z⃗)}(undef, length(μ∂²F_structure))
            Threads.@threads for t = 1:traj.T-1
                zₜ = Z⃗[slice(t, traj.dim)]
                zₜ₊₁ = Z⃗[slice(t + 1, traj.dim)]
                μₜ = μ⃗[slice(t, dynamics_dim)]
                μₜ∂²fₜ = μ∂²f(zₜ, zₜ₊₁, μₜ)
                for (k, (i, j)) ∈ enumerate(μ∂²f_structure)
                    μ∂²s[index(t, k, μ∂²f_nnz)] = μₜ∂²fₜ[i, j]
                end
            end
            return μ∂²s
        end
    else
        μ∂²F = nothing
    end

    return QuantumDynamics(
        integrators,
        F,
        ∂F,
        ∂F_structure,
        μ∂²F,
        μ∂²F_structure,
        dynamics_dim
    )
end

QuantumDynamics(P::AbstractIntegrator, traj::NamedTrajectory; kwargs...) =
    QuantumDynamics([P], traj; kwargs...)


end
