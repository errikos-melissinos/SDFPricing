import Base.Iterators as iter

xRanges = [1:2:3, 2:2:6]

myShape = (1, length.(xRanges)...)
means = Array{Vector}(undef, myShape)

j = 0
for indices in iter.product((length.(xRanges[end:-1:1]) .|> x -> 1:x)...)
    j += 1
    println(j)
    println(indices)
    means[j] = collect(indices)
end
means

for indices in iter.product(xRanges...)
    println(indices)
end


for indices in iter.product(xRanges[end:-1:1]...)
    println(indices)
end


u0 = vcat([[x, y, 0.0] for (x, y) in Base.Iterators.product(xRanges[1], xRanges[2])]...)

u0 = [(values..., 0.0) for values in Base.Iterators.product(xRanges...)]