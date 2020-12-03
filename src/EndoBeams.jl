module EndoBeams

# Include other files
include("packages.jl")
include("constants.jl")
include("rotations.jl")
include("shapefunctions.jl")
include("beamtypes.jl")
include("assembly.jl")
include("solver.jl")
include("contact.jl")

# Export the functions that you want to use outside of this packages
# Note: these functions should be commented using docstrings, which will help generating automatic documentation
export assemble!

end
