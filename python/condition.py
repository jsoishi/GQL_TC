import numpy as np
from mpi4py import MPI
import logging
logger = logging.getLogger(__name__)

def cond_number(solver):
    cond_nums = np.zeros(len(solver.pencils))
    for i,p in enumerate(solver.pencils):
        cond_nums[i] = np.linalg.cond(p.L.todense())
    comm = solver.domain.dist.comm
    global_max = np.array([0.,])
    local_max = np.array([cond_nums.max(),])
    comm.Reduce(local_max,global_max, MPI.MAX)
    if comm.rank == 0:
        logger.info("Maximum condition number is {:e}".format(global_max[0]))
