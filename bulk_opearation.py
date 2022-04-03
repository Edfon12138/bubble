import os
import time
def bulk(starter):
    for i in range(1,6):
        cmd="cd {}/{} && chmod +x runfile_feng && srun -p amd_512 -n 1 runfile_feng &".format(starter,i)
        os.system(cmd)
        print("done")
        time.sleep(2)

bulk("try")
