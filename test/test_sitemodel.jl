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


@testset "assign_rate_categories Function Tests" begin

    model = SiteModel(mutation_rate=0.1, gamma_shape=2.0, proportion_invariant=0.1, gamma_category_count=4, substitution_model=JC())

    # Test 1: Basic functionality
    @test length(assign_rate_categories(model, 100)) == 100

    # Test 2: Proportion invariant
    @test isapprox(sum(assign_rate_categories(model, 10_000) .== 1), 1_000; atol=500)

    # Test 3: Category range
    model_high_count = SiteModel(mutation_rate=0.1, gamma_shape=2.0, proportion_invariant=0.1, gamma_category_count=10, substitution_model=JC())
    categories = assign_rate_categories(model_high_count, 100)
    @test all(1 .<= categories .<= 11)  # 1 for invariant, 2:11 for variable

    # Test 6: Edge cases
    @test length(assign_rate_categories(model, 1)) == 1
    @test length(assign_rate_categories(model, 10000)) == 10000

    # Test 7: Input validation
    @test_throws ArgumentError assign_rate_categories(model, -100)

end
