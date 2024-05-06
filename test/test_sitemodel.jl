# First, test the proper creation and attribute setting of the SiteModel
@testset "SiteModel Construction Tests" begin
    substitution_model = JC()  # Assuming JC is a simple substitution model you have defined

    @testset "Valid construction" begin
        model = SiteModel(
            mutation_rate = 1.0,
            gamma_category_count = 4,
            gamma_shape = 0.5,
            proportion_invariant = 0.1,
            substitution_model = substitution_model
        )
        @test model.mutation_rate == 1.0
        @test model.gamma_category_count == 4
        @test model.gamma_shape == 0.5
        @test model.proportion_invariant == 0.1
        @test model.substitution_model === substitution_model
    end

    @testset "Invalid parameters" begin
        @test_throws ErrorException SiteModel(
            mutation_rate = -1.0,
            gamma_category_count = 4,
            gamma_shape = 0.5,
            proportion_invariant = 0.1,
            substitution_model = substitution_model
        )

        @test_throws ErrorException SiteModel(
            mutation_rate = 1.0,
            gamma_category_count = 0,  # Invalid category count
            gamma_shape = 0.5,
            proportion_invariant = 0.1,
            substitution_model = substitution_model
        )

        @test_throws ErrorException SiteModel(
            mutation_rate = 1.0,
            gamma_category_count = 4,
            gamma_shape = -0.5,  # Negative gamma shape
            proportion_invariant = 0.1,
            substitution_model = substitution_model
        )

        @test_throws ErrorException SiteModel(
            mutation_rate = 1.0,
            gamma_category_count = 4,
            gamma_shape = 0.5,
            proportion_invariant = 1.1,  # Invalid proportion invariant
            substitution_model = substitution_model
        )
    end
end
