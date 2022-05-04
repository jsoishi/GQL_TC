import numpy as np

def random_vector_potential(domain, R1, R2):

    r = domain.grid(-1, scales = domain.dealias)
    # Random perturbations to v in (r, z)
    gshape = domain.dist.grid_layout.global_shape(scales=domain.dealias)
    slices = domain.dist.grid_layout.slices(scales=domain.dealias)
    rand = np.random.RandomState(seed=42)

    filter_fraction = 0.5
    Ar = domain.new_field()
    Atheta = domain.new_field()
    Az = domain.new_field()

    fr= (4*(r-R1)*(R2-r))**4
    for A in [Ar, Atheta, Az]:
        A.set_scales(domain.dealias, keep_data=False)
        A['g'] = rand.standard_normal(gshape)[slices]
        A.set_scales(0.5,keep_data=True)
        A['c']
        A['g']
        A.set_scales(domain.dealias,keep_data=True)
        A['g'] *= fr

    return Ar, Atheta, Az
