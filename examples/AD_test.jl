using ChainRulesCore, Zygote
import ChainRulesCore: rrule
NO_FIELDS = ZeroTangent()

# Let's say we have z(y,x) = y(x)^2 , where y(x) = x+2

function f₁(x) # this will be z

    return sum(f₂(x).^2)

end

function f₂(x) # this will be y

    return x.^2.0.+5.0

end

function DyDx(x)

    return 2.0.*x

end

function rrule(::typeof(f₂), x)
    function y_pullback(dydx)
      NO_FIELDS, dydx.*DyDx(x)
    end
    f₂(x), y_pullback
end


x₀ = Float64[1.2, 3.4, 5.7]

grad = Zygote.gradient(x -> f₁(x), x₀)
println(grad[1])