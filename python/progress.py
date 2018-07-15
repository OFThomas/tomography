# Print progress to the console

import cProfile
import pstats
import numpy as np

def show_progress(pr,p):
    rem = "?"
    pr.disable()
    pr.create_stats()
    ps = pstats.Stats(pr)
    pr.enable()
    t = ps.total_tt
    if p != 0 : rem = str(np.around(t/p - t,3))
    print("\u001b[37mPercent done:\u001b[92m",str(np.around(100*p,2))+"%\u001b[37m",
          "\tTotal time:\u001b[94m",str(np.around(t,2))+"s\u001b[37m", 
          "\tEstimated time remaining:\u001b[95m",str(rem)+"s\u001b[37m")
    
