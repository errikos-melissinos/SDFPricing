import Base.Iterators as iter

xRanges = [1:2:3, 2:2:6]
myShape = (1, length.(xRanges)...)
means = Array{Vector}(undef, myShape)

j = 0
for indices in iter.product((length.(xRanges[end:-1:1]) .|> x -> 1:x)...)
    j += 1
    println(j)
    println(indices)
    means = collect(indices)
end
