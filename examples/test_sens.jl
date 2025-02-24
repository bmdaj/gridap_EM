using Zygote

function FOM_Ez(Ez_vec)


    val = Ez_vec' * Ez_vec
    return real(val)

end

N = 100000

Ez_vec = rand(ComplexF64, N)

grad = Zygote.gradient(x -> FOM_Ez(x), Ez_vec)

grad[1]