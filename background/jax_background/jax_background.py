from cosmosis.datablock import option_section, names
import jax_cosmo
import jax.numpy as jnp

cosmosis_jax = True

def setup(options):
    zmax = options.get_double(option_section, "zmax", default=3.0)
    nz = options.get_int(option_section, "nz", default=301)
    return {
        "zmax": zmax,
        "nz": nz,
    }

inputs = [
    (names.cosmological_parameters, "omega_m"),
    (names.cosmological_parameters, "omega_k"),
    (names.cosmological_parameters, "w"),
    (names.cosmological_parameters, "h0"),
]

def execute(inputs, config):
    omega_m = inputs[("cosmological_parameters", "omega_m")]
    omega_k = inputs[("cosmological_parameters", "omega_k")]
    w = inputs[("cosmological_parameters", "w")]
    h0 = inputs[("cosmological_parameters", "h0")]
    omega_b = 0.044
    omega_c = omega_m - omega_b
    cosmo = jax_cosmo.Cosmology(Omega_c=omega_c, Omega_b=omega_b, h=h0, n_s=0.96, sigma8=0.8, Omega_k=omega_k, w0=w, wa=0.0)
    z = jnp.linspace(0.0, config["zmax"], config["nz"])[1:]
    a = 1 / (1 + z)
    d_a = jax_cosmo.background.angular_diameter_distance(cosmo, a) / h0
    d_l = (1 + z)**2 * d_a
    mu = 5 * jnp.log10(d_l) + 25
    return {
        ("distances", "z"): z,
        ("distances", "a"): a,
        ("distances", "mu"): mu,
        ("distances", "d_a"): d_a,
        ("distances", "d_l"): d_l,
    }
