module Integrators

export AbstractIntegrator

export QuantumIntegrator

export QuantumPadeIntegrator
export QuantumStatePadeIntegrator
export UnitaryPadeIntegrator

export UnitaryExponentialIntegrator
export QuantumStateExponentialIntegrator
export DensityOperatorExponentialIntegrator

export DerivativeIntegrator

export jacobian
export hessian_of_the_lagrangian

export get_comps

export nth_order_pade
export fourth_order_pade
export sixth_order_pade
export eighth_order_pade
export tenth_order_pade

using NamedTrajectories
using TrajectoryIndexingUtils
using PiccoloQuantumObjects
using LinearAlgebra
using SparseArrays
using ForwardDiff
using TestItems

const ⊗ = kron

abstract type AbstractIntegrator end

function comps(P::AbstractIntegrator, traj::NamedTrajectory)
    state_comps = traj.components[state(P)]
    u = controls(P)
    if u isa Tuple
        control_comps = Tuple(traj.components[uᵢ] for uᵢ ∈ u)
    else
        control_comps = traj.components[u]
    end
    if traj.timestep isa Float64
        return state_comps, control_comps
    else
        timestep_comp = traj.components[traj.timestep]
        return state_comps, control_comps, timestep_comp
    end
end



abstract type QuantumIntegrator <: AbstractIntegrator end

abstract type UnitaryIntegrator <: QuantumIntegrator end
abstract type QuantumStateIntegrator <: QuantumIntegrator end
abstract type DensityOperatorIntegrator <: QuantumIntegrator end

function update_state_components!(
    integrator::QuantumIntegrator, 
    state_name::Symbol,
    traj::NamedTrajectory
)
    integrator.state_components = traj.components[state_name]
    return nothing 
end

function update_drive_components!(
    integrator::QuantumIntegrator, 
    drive_name::Symbol,
    traj::NamedTrajectory
)
    integrator.drive_components = traj.components[drive_name]
    return nothing 
end

function update_variable_components!(
    integrator::DerivativeIntegrator,
    variable_name::Symbol,
    derivative_name::Symbol,
    traj::NamedTrajectory
)
    integrator.variable_components = traj.components[variable_name]
    return nothing 
end



include("_integrator_utils.jl")
include("derivative_integrator.jl")
include("pade_integrators.jl")
include("exponential_integrators.jl")


end
