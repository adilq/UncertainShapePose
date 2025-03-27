import numpy as np
import matplotlib.pyplot as plt
from itertools import cycle

step_data = np.loadtxt(open("Steps_dynamic_example.csv", "rb"), delimiter=",")
print(step_data)

num_obj = 10
num_step = 50
stationarity = np.zeros((num_obj,num_step))
for i in range(step_data.shape[0]):
    row = step_data[i,:]
    step = int(row[0])
    id = int(row[1])
    score = row[2]
    gt_state = int(row[3])
    stationarity[id,step] = score

labels = ['Obj 0','Obj 1','Obj 2','Obj 3','Obj 4','Obj 5','Obj 6','Obj 7','Obj 8','Obj 9']
styles = ['^--', '-' , '-', '-', '-', 's--', 'p--', '*--', '.--', 'v--']
linecycler = cycle(styles)

steps = [i for i in range(num_step+1)]
plt.figure()
for i in range(stationarity.shape[0]):
    conf = stationarity[i,:]
    plt.plot(steps, np.concatenate(([0.5],conf)), next(linecycler)[::-1], linewidth=0.3, label=labels[i])
plt.xlim((0,10))
plt.ylim((0.2,0.85))
plt.xlabel('Step')
plt.ylabel('E[v]')
plt.legend()
plt.title("Object Level Confidence E[v] Over Steps")
#plt.show()
plt.savefig('steps.png', dpi=600)
