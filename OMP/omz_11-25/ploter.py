import numpy as np
import matplotlib.pyplot as plt
from numpy import arange


Thread = [2,4,6,8,10,12]
Told125 = 951.861539
Tnew125 = [584.994349, 327.380138, 260.283109, 250.63611, 255.046343, 227.633683]

Told251 = 7652.136362
Tnew251 = [2164.829393, 1821.036122, 1551.989532, 1458.806854, 1701.936160, 1436.028891]

speedup125 = []
speedup251 = []
efficiency125 = []
efficiency251 = []

for i in range(6):
    speedup125.append(Told125/Tnew125[i])
    efficiency125.append(speedup125[i]/Thread[i])

for i in range(6):
    speedup251.append(Told251/Tnew251[i])
    efficiency251.append(speedup251[i]/Thread[i])

plt.plot(Thread, speedup125, 'blue', label='125x125x125')
plt.plot(Thread, speedup251, 'red', label='251x251x251')


'''fig, ax = plt.subplots()
fig.canvas.draw()
ai_list=[1/16., 1/8., 1/4., 1/2., 1, 2, 3, 4, 6, 8, 16]'''
plt.legend( loc='upper right')
plt.ylabel('Speedup')
plt.xlabel('Thread')
plt.title('Speedup')
plt.grid(True)
plt.savefig('speedup_3Domp.png')
plt.show()

plt.plot(Thread, efficiency125, 'blue', label='125x125x125')
plt.plot(Thread, efficiency251, 'red', label='251x251x251')


'''fig, ax = plt.subplots()
    fig.canvas.draw()
ai_list=[1/16., 1/8., 1/4., 1/2., 1, 2, 3, 4, 6, 8, 16]'''
plt.legend( loc='upper right')
plt.ylabel('Efficiency')
plt.xlabel('Thread')
plt.title('Efficiency')
plt.grid(True)
plt.savefig('efficiency_3Domp.png')
plt.show()





