using UCIWWEIHR
using Test

@testset "UCIWWEIHR.jl" begin
    # Write your tests here.
    x = 2
    y = 2
    @test UCIWWEIHR.sum_values(x,y) == 4
    
end
