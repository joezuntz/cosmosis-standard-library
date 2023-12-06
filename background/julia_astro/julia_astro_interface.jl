module julia_astro_interface
src_dir = ENV["COSMOSIS_SRC_DIR"]
julia_dir = "$(src_dir)/datablock/julia"
println(julia_dir)
push!(LOAD_PATH, julia_dir)
using cosmosis
import Cosmology
import ForwardDiff



function setup(options::cosmosis.datablock)

    zmax = cosmosis.get_double(options, cosmosis.option_section, "zmax")
    nz = cosmosis.get_int(options, cosmosis.option_section, "nz")
    d = Dict()
    d["zmax"] = zmax
    d["nz"] = nz
    return d

end

inputs = [
    ("cosmological_parameters", "omega_m"),
    ("cosmological_parameters", "h0"),
    ("cosmological_parameters", "w"),
    ("cosmological_parameters", "omega_k")
]

# Dict{Tuple{String,String},   Union{T, Array{T}}  }

function execute(inputs::Dict{Tuple{String,String}, Union{T, Array{T}}}, config::Dict) where T
    # Extract the inputs we need from the relevant section
    section = "cosmological_parameters"
    omega_m = inputs[section, "omega_m"]
    h0 = inputs[section, "h0"]
    w = inputs[section, "w"]
    omega_k = inputs[section, "omega_k"] 

    # Define the z samples we want
    zmax = config["zmax"]
    nz = config["nz"]
    dz = zmax / (nz - 1)

    # Redshift range
    z = LinRange(0.0, zmax, nz)

    # Main calculation
    # We cannot use the Cosmology.cosmology function because it chooses the
    # cosmology based on the parameters, and if you pick e.g. omega_k=0 or w=-1
    # then it generates a subclass that ignores those parameters. Then the
    # derivative calculation doesn't work on those parameters because they are
    # ignored.
    omega_lambda = 1 - omega_k - omega_m
    omega_r = 0.0
    wa = 0.0
    c = Cosmology.WCDM(h0, omega_k, omega_lambda, omega_m, omega_r, w, wa)

    # The Cosmology module redshift input does not work for a generic type T
    # array, because of the way the underlying integration is implemented. It
    # also looks too painful to fix.  But the z is fixed anyway, so we can just
    # make a normal Float64 array of z values here and use those to calculate
    # the distances.  Then we can replace the z with the T version. It won't have
    # derivative information in, but that's okay as it's constant anyway
    z = collect(0.0:dz:zmax)
    d_a = [Cosmology.angular_diameter_dist(c, zi).val for zi in z]
    mu = [Cosmology.distmod(c, zi) for zi in z]

    # Putting the actual first value of -inf in the array causes problems for
    # ForwardDiff, so we just put a zero in there instead
    mu[1] = T(0.0)

    # now replace with the T version, just to make ForwardDiff happy
    z = Vector{T}(z)

    return Dict(
        ("distances", "z") => z,
        ("distances", "d_a") => d_a,
        ("distances", "mu") => mu,
    )

end



end
    