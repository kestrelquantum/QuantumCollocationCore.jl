
###
### Derivative Integrator
###

mutable struct DerivativeIntegrator <: AbstractIntegrator
    variable_components::Vector{Int}
    derivative_components::Vector{Int}
    freetime::Bool
    timestep::Union{Float64, Int} # timestep of index in z
    zdim::Int
    dim::Int
    autodiff::Bool
end

function DerivativeIntegrator(
    variable::Symbol,
    derivative::Symbol,
    traj::NamedTrajectory
)
    freetime = traj.timestep isa Symbol
    if freetime
        timestep = traj.components[traj.timestep][1]
    else
        timestep = traj.timestep
    end
    return DerivativeIntegrator(
        traj.components[variable],
        traj.components[derivative],
        freetime,
        timestep,
        traj.dim,
        traj.dims[variable],
        false
    )
end

# function (integrator::DerivativeIntegrator)(
#     traj::NamedTrajectory;
#     variable::Union{Symbol, Nothing}=nothing,
#     derivative::Union{Symbol, Nothing}=nothing,
# )
#     @assert !isnothing(variable) "variable must be provided"
#     @assert !isnothing(derivative) "derivative must be provided"
#     return DerivativeIntegrator(
#         variable,
#         derivative,
#         traj
#     )
# end

state(integrator::DerivativeIntegrator) = integrator.variable
controls(integrator::DerivativeIntegrator) = integrator.derivative

@views function (D::DerivativeIntegrator)(
    zₜ::AbstractVector,
    zₜ₊₁::AbstractVector,
    t::Int
)
    xₜ = zₜ[D.variable_components]
    xₜ₊₁ = zₜ₊₁[D.variable_components]
    dxₜ = zₜ[D.derivative_components]
    if D.freetime
        Δtₜ = zₜ[D.timestep]
    else
        Δtₜ = D.timestep
    end
    return xₜ₊₁ - xₜ - Δtₜ * dxₜ
end

@views function jacobian(
    D::DerivativeIntegrator,
    zₜ::AbstractVector,
    zₜ₊₁::AbstractVector,
    t::Int
)
    ∂D = spzeros(D.dim, 2D.zdim)
    if D.freetime
        dxₜ = zₜ[D.derivative_components]
        Δtₜ = zₜ[D.timestep]
        ∂D[:, D.timestep] .= -dxₜ
    else
        Δtₜ = D.timestep
    end
    ∂D[:, D.variable_components] .= -1.0I(D.dim)
    ∂D[:, D.derivative_components] .= -Δtₜ * I(D.dim)
    ∂D[:, D.zdim .+ D.variable_components] .= 1.0I(D.dim)
    return ∂D
end

function get_comps(D::DerivativeIntegrator, traj::NamedTrajectory)
    if D.freetime
        return D.variable_components, D.derivative_components, traj.components[traj.timestep]
    else
        return D.variable_components, D.derivative_components
    end
end
