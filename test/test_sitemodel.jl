@testset "site model invariants" begin
    Random.seed!(11)
    site_model = SiteModel(100, 0.01, 4, 0.5, 0.2, JC())

    @test site_model.sequence_length == 100
    @test length(site_model.variable_sites) == 4
    @test length(site_model.μ) == 4
    @test all(rate -> isfinite(rate) && rate > 0.0, site_model.μ)

    variable_sites = reduce(vcat, site_model.variable_sites)
    @test length(variable_sites) == 80
    @test sort(variable_sites) == unique(sort(variable_sites))
    @test all(1 .<= variable_sites .<= 100)
    @test isapprox(length.(site_model.variable_sites) ⋅ site_model.μ / site_model.sequence_length, site_model.mutation_rate; rtol=1e-10)
end

@testset "site model without gamma rate categories" begin
    Random.seed!(12)
    site_model = SiteModel(10, 0.5, 0, 0.0, 0.3, JC())
    @test length(site_model.variable_sites) == 1
    @test length(site_model.variable_sites[1]) == 7
    @test site_model.μ == [0.5]
end

@testset "site model validation" begin
    @test_throws ArgumentError SiteModel(0, 0.1, 1, 1.0, 0.0, JC())
    @test_throws ArgumentError SiteModel(10, 0.0, 1, 1.0, 0.0, JC())
    @test_throws ArgumentError SiteModel(10, 0.1, -1, 1.0, 0.0, JC())
    @test_throws ArgumentError SiteModel(10, 0.1, 1, 0.0, 0.0, JC())
    @test_throws ArgumentError SiteModel(10, 0.1, 1, 1.0, -0.1, JC())
    @test_throws ArgumentError SiteModel(10, 0.1, 1, 1.0, 1.0, JC())
    @test_throws ArgumentError SeqSim.assign_rates(10, 0.0, 0.1, 0.0, 1)
end
