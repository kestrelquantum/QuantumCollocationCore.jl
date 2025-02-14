module Integrators

export AbstractIntegrator

export QuantumIntegrator

export jacobian
export hessian_of_the_lagrangian

export get_comps

using NamedTrajectories
using TrajectoryIndexingUtils
using PiccoloQuantumObjects
using LinearAlgebra
using SparseArrays
using ForwardDiff
using TestItems

import .NamedTrajectories: add_suffix, remove_suffix, get_suffix

const ⊗ = kron

"""
    AbstractIntegrator

Abstract type for integrators.
    
Required methods:
    - Integrator(args): evaluate an implicit dynamics constraint between knot points
    - jacobian(Integrator, args): evaluate the jacobian of the implicit dynamics

Optional methods:
    - hessian_of_the_lagrangian: evaluate the hessian of the lagrangian

"""
abstract type AbstractIntegrator end

function jacobian end

function hessian_of_the_lagrangian end

abstract type QuantumIntegrator <: AbstractIntegrator end

abstract type UnitaryIntegrator <: QuantumIntegrator end
abstract type QuantumStateIntegrator <: QuantumIntegrator end
abstract type DensityOperatorIntegrator <: QuantumIntegrator end

include("_integrator_utils.jl")
include("derivative_integrator.jl")
include("pade_integrators.jl")
include("exponential_integrators.jl")


# ---------------------------------------------------------------------------- #
# Integrator add/remove/get suffix methods
# ---------------------------------------------------------------------------- #

# TODO: Refactor (don't overload, use get_comps)


"""
    modify_integrator_suffix(integrator, modifier, suffix, traj, mod_traj)

Modify the integrator by adding or removing a suffix to the component names. 
The `modifier` function should be either `add_suffix` or `remove_suffix`.

Integrators use components and not symbols, so the suffix is added or removed by using 
the components of modifier(name, suffix) from `mod_traj`.

"""
function modify_integrator_suffix(
    integrator::AbstractIntegrator,
    modifier::Function,
    suffix::String,
    traj::NamedTrajectory,
    mod_traj::NamedTrajectory
)
    mod_integrator = deepcopy(integrator)

    if integrator isa QuantumIntegrator
        state_name = get_component_names(traj, integrator.state_components)
        drive_name = get_component_names(traj, integrator.drive_components)
        mod_integrator.state_components = mod_traj.components[modifier(state_name, suffix)]
        mod_integrator.drive_components = mod_traj.components[modifier(drive_name, suffix)]
        return mod_integrator
    elseif integrator isa DerivativeIntegrator
        var_name = get_component_names(traj, integrator.variable_components)
        der_name = get_component_names(traj, integrator.derivative_components)
        mod_integrator.variable_components = mod_traj.components[modifier(var_name, suffix)]
        mod_integrator.derivative_components = mod_traj.components[modifier(der_name, suffix)]
        return mod_integrator
    else
        error("Integrator type not recognized")
    end
end

# add suffix

function add_suffix(
    integrator::AbstractIntegrator,
    suffix::String,
    traj::NamedTrajectory,
    mod_traj::NamedTrajectory
)
    return modify_integrator_suffix(integrator, add_suffix, suffix, traj, mod_traj)
end

function add_suffix(
    integrators::AbstractVector{<:AbstractIntegrator},
    suffix::String,
    traj::NamedTrajectory,
    mod_traj::NamedTrajectory
)
    return [add_suffix(i, suffix, traj, mod_traj) for i ∈ integrators]
end

# remove suffix

function remove_suffix(
    integrator::AbstractIntegrator,
    suffix::String,
    traj::NamedTrajectory,
    mod_traj::NamedTrajectory
)
    return modify_integrator_suffix(integrator, remove_suffix, suffix, traj, mod_traj)
end

function remove_suffix(
    integrators::AbstractVector{<:AbstractIntegrator},
    suffix::String,
    traj::NamedTrajectory,
    mod_traj::NamedTrajectory
)
    return [remove_suffix(i, suffix, traj, mod_traj) for i in integrators]
end


# get suffix

function get_suffix(
    integrators::AbstractVector{<:AbstractIntegrator},
    traj::NamedTrajectory,
    mod_traj::NamedTrajectory,
    suffix::String;
    remove::Bool=false
)
    found = AbstractIntegrator[]
    for intg ∈ integrators
        if intg isa QuantumIntegrator
            name = get_component_names(traj, intg.state_components)
        elseif intg isa DerivativeIntegrator
            name = get_component_names(traj, intg.variable_components)
        else
            error("Integrator type not recognized")
        end

        if endswith(name, suffix)
            if remove
                push!(found, remove_suffix(intg, suffix, traj, mod_traj))
            else
                push!(found, intg)
            end
        end
    end
    return found
end


end
