# Print progress to the console

import cProfile
import pstats

def show_progress(pr,p):
    rem = "?"
    pr.disable()
    pr.create_stats()
    ps = pstats.Stats(pr)
    pr.enable()
    t = ps.total_tt
    if p != 0 : rem = str(round(t/p - t,3))
    print("\u001b[0m\u001b[1mPercent done:\u001b[32m",str(round(100*p,2))+"%\u001b[0m\u001b[1m",
          "\tTotal time:\u001b[34m",str(round(t,2))+"s\u001b[0m\u001b[1m", 
          "\tEstimated time remaining:\u001b[35m",str(rem)+"s\u001b[0m\u001b[1m",
          end='                 \r')
    
