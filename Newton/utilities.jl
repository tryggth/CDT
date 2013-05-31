# CDT Utilities
#
# Various utility functions ported from cdt-utilities.lisp
# List operations in Lisp are replaced by tuple or 1-d array operations
# in Julia.
# Only dimension independent utilities are in this file
# Many functions such as standard deviation and mean are built into Julia already,
# and hence do not need to be defined here.

function sum(list)
	# Sums the elements in a tuple or 1-d array
	return apply(+, list)
end

