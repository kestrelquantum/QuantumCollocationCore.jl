export DensityOperatorPureStateInfidelityLoss

struct DensityOperatorPureStateInfidelityLoss <: AbstractLoss
    l::Function
    ∇l::Function
    ∇²l::Function
    ∇²l_structure::Vector{Tuple{Int,Int}}
    name::Symbol

    function DensityOperatorPureStateInfidelityLoss(
        name::Symbol,
        ψ_goal::AbstractVector
    )
        ψ_real = real(ψ_goal)
        ψ_imag = imag(ψ_goal)
        ψ_goal_lifted = vcat(
            ψ_real ⊗ ψ_real + ψ_imag ⊗ ψ_imag,
            ψ_real ⊗ ψ_imag - ψ_imag ⊗ ψ_real
        )

        l = ρ⃗̃ -> abs(1 - ψ_goal_lifted' * ρ⃗̃)

        ∇l = ρ⃗̃ -> -sign(1 - ψ_goal_lifted' * ρ⃗̃) * ψ_goal_lifted

        ∇²l = ρ⃗̃ -> spzeros(length(ρ⃗̃), length(ρ⃗̃))
        ∇²l_structure = []

        return new(l, ∇l, ∇²l, ∇²l_structure, name)
    end
end

function (loss::DensityOperatorPureStateInfidelityLoss)(
    ρ⃗̃::AbstractVector{<:Real};
    gradient=false,
    hessian=false
)
    @assert !(gradient && hessian)
    if !(gradient || hessian)
        return loss.l(ρ⃗̃)
    elseif gradient
        return loss.∇l(ρ⃗̃)
    elseif hessian
        return loss.∇²l(ρ⃗̃)
    end
end

# TODO: write tests
