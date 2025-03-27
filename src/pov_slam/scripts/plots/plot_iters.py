import numpy as np
import matplotlib.pyplot as plt
from matplotlib.ticker import MaxNLocator
from itertools import cycle

step_data = np.loadtxt(open("Iterations_dynamic_example.csv", "rb"), delimiter=",")
print(step_data)

num_obj = 10
num_step = 50
num_iter = 100
stationarity = np.zeros((num_obj,num_step,num_iter))
rin = np.zeros((num_obj,num_step,num_iter))
for i in range(step_data.shape[0]):
    row = step_data[i,:]
    step = int(row[0])
    iter = int(row[1])
    id = int(row[2])
    score = row[3]
    r_score = row[4]
    gt_state = int(row[5])
    stationarity[id,step,iter] = score
    rin[id,step,iter] = r_score

labels = ['Obj 0','Obj 1','Obj 2','Obj 3','Obj 4','Obj 5','Obj 6','Obj 7','Obj 8','Obj 9']
styles = ['^--', '-' , '-', '-', '-', 's--', 'p--', '*--', '.--', 'v--']
linecycler = cycle(styles)

num_iter = 20
target_step = 0
iters = [int(i) for i in range(num_iter+1)]
plt.figure()
for i in range(stationarity.shape[0]):
    conf = stationarity[i,target_step,0:num_iter]
    plt.plot(iters,np.concatenate(([0.5],conf)), next(linecycler)[::-1], linewidth=0.6, label=labels[i])
plt.xlim((0,21))
plt.ylim((0.0,1))
plt.legend()
plt.title("Object Level Confidence Over Iteration at Frame "+str(target_step))
#plt.show()
plt.xticks(np.arange(0, num_iter+1, step=4))
plt.xlabel('Iteration')
plt.ylabel('Stationarity')
plt.savefig('stationarity.png', dpi=300)


plt.figure()
for i in range(rin.shape[0]):
    conf = rin[i,target_step,0:num_iter]
    plt.plot(iters,np.concatenate(([0.5],conf)), next(linecycler)[::-1], linewidth=0.6, label=labels[i])
plt.xlim((0,21))
plt.ylim((0.0,1))
plt.legend()
plt.title("Object Level E[pi] Over Iteration at Frame "+str(target_step))
#plt.show()
plt.xticks(np.arange(0, num_iter+1, step=4))
plt.xlabel('Iteration')
plt.ylabel('E[pi]')
plt.savefig('rin.png', dpi=300)


errors = np.loadtxt(open("pose_error.csv", "rb"), delimiter=",")
frame_0 = errors[0:num_iter+1,0]

plt.figure()
plt.plot(iters,frame_0, linewidth=0.6)
plt.xlim((0,num_iter+1))
plt.ylim((0.0,0.28))
plt.title("Robot Pose Error Over Iteration at Frame 0")
plt.xticks(np.arange(0, num_iter+1, step=4))
plt.xlabel('Iteration')
plt.ylabel('Error')
plt.savefig('pose_error.png', dpi=300)
