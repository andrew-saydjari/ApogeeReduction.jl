using SparseArrays, LinearAlgebra

function I_twiddle(n::Int)
    view(I(n), n:-1:1, :)
end

function construct_P(n)
    P = spzeros(Bool, n, n)
    for i in 1:2:n
        P[(i + 1) ÷ 2, i] = true
    end
    for i in 2:2:n
        P[(n + 1) ÷ 2 + (i ÷ 2), i] = true
    end
    P
end

function nsin2(n::Int)
    if n == 2
        [1.0 1.0
         1.0 -1.0]
    elseif n >= 4
        n1 = n ÷ 2
        M1 = nsin4(n1)
        M2 = nsin2(n1)

        I_n1 = I(n1)
        Ĩ = I_twiddle(n1)

        H = (1 / sqrt(2)) .* [I(n1) Ĩ
                              I(n1) -Ĩ]
        P = construct_P(n)
        P' * [M1 zeros(n1, n1)
              zeros(n1, n1) M2] * 2H
    end
end

function nsin3(n::Int)
    if n == 2
        return [1.0 1.0; 1.0 -1.0]
    elseif n >= 4
        n1 = n ÷ 2
        M1 = nsin4(n1)
        M2 = nsin3(n1)
        I_n1 = I(n1)
        sqrt2H_transpose = [I_n1 I_n1; view(I_n1, n1:-1:1, :) -view(I_n1, n1:-1:1, :)]
        P = [I_n1 zeros(n1, n1); zeros(n1, n1) view(I_n1, n1:-1:1, :)]
        return sqrt2H_transpose * [M1 zeros(n1, n1); zeros(n1, n1) M2] * P
    end
end

#this is broken, but it's not needed for nsin1
function nsin4(n::Int)
    if n == 2
        sqrt(2) * [sin(π / 8) cos(π / 8)
                   cos(π / 8) -sin(π / 8)]
    elseif n >= 4
        n1 = n ÷ 2
        M1 = nsin2(n1) # TODO this is not a typo??? It's like this in the paper.
        M2 = nsin2(n1)
        I_n1 = I(n1)
        V = [zeros(1, n1 - 1) sqrt(2)
             I_n1[1:(n1 - 1), :] -I_n1[1:(n1 - 1), :]
             -I_n1[1:(n1 - 1), :] -I_n1[1:(n1 - 1), :]
             zeros(1, n1 - 1) sqrt(2)] / sqrt(2)
        @show V
        Q = [diagm(0 => [sin((2k + 1)π / (4n)) for k in 0:(n1 - 1)]) zeros(n1, n1);
             zeros(n1, n1) diagm(0 => [-cos((2k + 1)π / (4n)) for k in 0:(n1 - 1)])] *
            [I_n1 view(I_n1, n1:-1:1, :); -view(I_n1, n1:-1:1, :) I_n1]
        @show Q
        P = [I_n1 zeros(n1, n1); zeros(n1, n1) view(I_n1, n1:-1:1, :)]
        @show P
        P' * (sqrt(2) * V) * [M1 zeros(n1, n1); zeros(n1, n1) M2] * Q
    end
end

"""
This is a function n, although it's written as nsin1(n-1) in the paper
"""
function nsin1(n::Int)
    #assert that n is a power of 2
    if ((n & (n - 1)) != 0) || (n <= 0)
        throw(ArgumentError("n must be a power of 2"))
    end
    if n == 2
        return [1.0]
    elseif n >= 4
        n1 = n ÷ 2 # this definition from paper

        # construct H
        Ĩ = I_twiddle(n1 - 1)
        H_hat = [I(n1 - 1) zeros(n1 - 1, 1) Ĩ;
                 zeros(1, n1 - 1) sqrt(2) zeros(1, n1 - 1);
                 I(n1 - 1) zeros(n1 - 1, 1) -Ĩ] / sqrt(2)

        P = construct_P(n - 1)
        M1 = nsin3(n1)
        M2 = nsin1(n1)

        return P' * [M1 zeros(n1, n1 - 1); zeros(n1 - 1, n1) M2] * (sqrt(2) * H_hat)
    end
end