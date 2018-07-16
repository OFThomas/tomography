# Contains a function for printing progress

import cProfile
import pstats
import numpy as np

# p = proportion done
# pr is a profiling object
def show_progress(pr, p):
    rem = "?"
    pr.disable()
    pr.create_stats()
    ps = pstats.Stats(pr)
    t = ps.total_tt
    pr.enable()
    if(p != 0) : rem = str(np.around(t/p-t,2))
    print("Percent done:\u001b[93m", str(np.around(100*p,3))+"%\u001b[37m",
          "\tTotal time:\u001b[91m",str(np.around(t,1))+"s\u001b[37m",
          "\tEstimated time remaining:\u001b[92m",rem+"s\u001b[37m")
