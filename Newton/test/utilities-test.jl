# Test of CDT Utilities
#
# This provides various tests for the utility functions in utilities.jl based on
# the [FactCheck] framework. To install do 'Pkg.add("FactCheck") at the julia prompt

require("../utilities.jl")

using FactCheck
# Allows functional testing

@facts "Utility functions work" begin

	x = (1,2,3,4,5,6)
	y = [1,2,3,4,5,6]
	
	@fact "sum(list) works on tuples and 1-d arrays in" begin 
	
		sum(x) => 21
		# Tuples sum correctly
	
		sum(y) => 21
		# 1-d arrays sum correctly
	end
	
	@fact "Mean is already defined in Julia" begin
	
		mean(x)=>3.5
		# Tuples
		
		mean(y)=>3.5
		# Arrays
	end
	
	@fact "Standard deviation is already defined in Julia" begin
	
		std(y)=>1.8708286933869707
		# Arrays or vectors only
		
	end
		
end

# ## References
#
# [FactCheck]: https://github.com/zachallaun/FactCheck.jl