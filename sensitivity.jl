using ChainRulesCore, Zygote
import ChainRulesCore: rrule
NO_FIELDS = ZeroTangent()

dξt_dξf(ξ_f, β, η) = β * (1.0 - tanh(β * (ξ_f - η))^2) / (tanh(β * η) + tanh(β * (1.0 - η)))

dε_dξt(dξ_f, β, η) = (ε1 - ε0) * dξt_dξf(ξ_f, β, η)

