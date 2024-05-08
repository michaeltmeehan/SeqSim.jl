# First, test the proper creation and attribute setting of the SiteModel
@testset "SiteModel Construction Tests" begin
    substitution_model = JC()  # Assuming JC is a simple substitution model you have defined

    @testset "Default parameters" begin
        model = SiteModel(1e-3)
        @test model.mutation_rate == 1e-3
        @test model.gamma_shape == 0.0
        @test model.proportion_invariant == 0.0
        @test typeof(model.substitution_model) == JC
    end

    @testset "Valid construction" begin
        model = SiteModel(
            mutation_rate = 1.0,
            gamma_shape = 0.5,
            proportion_invariant = 0.1,
            substitution_model = substitution_model
        )
        @test model.mutation_rate == 1.0
        @test model.gamma_shape == 0.5
        @test model.proportion_invariant == 0.1
        @test model.substitution_model === substitution_model
    end

    @testset "Invalid parameters" begin
        @test_throws ErrorException SiteModel(
            mutation_rate = -1.0,
            gamma_shape = 0.5,
            proportion_invariant = 0.1,
            substitution_model = substitution_model
        )

        @test_throws ErrorException SiteModel(
            mutation_rate = 1.0,
            gamma_shape = -0.5,  # Negative gamma shape
            proportion_invariant = 0.1,
            substitution_model = substitution_model
        )

        @test_throws ErrorException SiteModel(
            mutation_rate = 1.0,
            gamma_shape = 0.5,
            proportion_invariant = 1.1,  # Invalid proportion invariant
            substitution_model = substitution_model
        )
    end
end



# Define a test suite for site_rates
@testset "SiteRates Function Tests" begin

    # Test 1: Basic functionality
    model = SiteModel(mutation_rate=0.1, gamma_shape=2.0, proportion_invariant=0.1, substitution_model=JC())
    rates = site_rates(model, 100_000)
    @test length(rates) == 100_000
    @test isapprox(sum(rates .== 0.0), 10_000; atol=500)  # Allowing some variation
    @test isapprox(mean(rates), 0.1; atol=1e-3)

    # Test 2: Zero gamma shape
    model_zero_gamma = SiteModel(mutation_rate=0.1, gamma_shape=0.0, proportion_invariant=0.1, substitution_model=JC())
    rates_zero_gamma = site_rates(model_zero_gamma, 100)
    @test all(rates_zero_gamma .== 0.1)

    # Test 3: Edge cases for proportion_invariant
    model_no_inv = SiteModel(mutation_rate=0.1, gamma_shape=2.0, proportion_invariant=0.0, substitution_model=JC())
    rates_no_inv = site_rates(model_no_inv, 100)
    @test all(rates_no_inv .> 0.0)

    # Test 4: Input validation
    @test_throws ArgumentError site_rates(model, -1)

end

