export DensityOperatorPureStateInfidelityObjective

function DensityOperatorPureStateInfidelityObjective(
    name::Symbol,
    goal::AbstractVector{<:Real};
    Q::Float64=100.0,
    eval_hessian::Bool=true
)
    loss = :DensityOperatorPureStateInfidelityLoss
    l = eval(loss)(name, goal)

    params = Dict(
        :type => loss,
        :name => name,
        :goal => goal,
        :Q => Q,
        :eval_hessian => eval_hessian
    )

    @views function L(Z⃗::AbstractVector{<:Real}, Z::NamedTrajectory)
        return Q * l(Z⃗[slice(Z.T, Z.components[name], Z.dim)])
    end

    @views function ∇L(Z⃗::AbstractVector{<:Real}, Z::NamedTrajectory)
        ∇ = zeros(Z.dim * Z.T)
        ρ⃗̃_slice = slice(Z.T, Z.components[name], Z.dim)
        ρ⃗̃ = Z⃗[ρ⃗̃_slice]
        ∇l = l(ρ⃗̃; gradient=true)
        ∇[ρ⃗̃_slice] = Q * ∇l
        return ∇
    end

    ∂²L_structure(Z::NamedTrajectory) = []
    ∂²L(Z⃗::AbstractVector{<:Real}, Z::NamedTrajectory) = []

	return Objective(L, ∇L, ∂²L, ∂²L_structure, Dict[params])
end
