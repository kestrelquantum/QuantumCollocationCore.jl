module Losses

export Loss

using NamedTrajectories
using TrajectoryIndexingUtils
using PiccoloQuantumObjects
using LinearAlgebra
using SparseArrays
using ForwardDiff
using TestItems
using ExponentialAction

# TODO:
# - [ ] Do not reference the Z object in the loss (components only / remove "name")

# ----------------------------------------------------------------------------- #
#                           Abstract Loss                                       #
# ----------------------------------------------------------------------------- #

abstract type AbstractLoss end

include("_experimental_loss_functions.jl")
include("quantum_state_infidelity_loss.jl")
include("unitary_trace_loss.jl")
include("unitary_infidelity_loss.jl")
include("density_operator_losses.jl")

end
