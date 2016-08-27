import numpy as np
from LDPC import *

H=LDPC('gallager_regular')
res=H.simBEC(0.4,50,10000000,100,1000,1000)

print '''Iterations %d BERs: %d FERs: %d BER%% : %f FER%% : %f'''%res
