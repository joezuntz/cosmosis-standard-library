import numpy as np
import matplotlib.pyplot as plt 



fig=plt.figure(figsize=(10,10))
ax=fig.add_subplot(111)
ax.set_xlim(-4,1)
ax.set_ylim(-4,1)

fig.patch.set_visible(False)
ax.axis('off')




ax.arrow(-3.5, 0, 4,0,lw=3, head_width=0.08, head_length=0.1, fc='k', ec='k')
ax.arrow(0, -3.5, 0,4,lw=3, head_width=0.08, head_length=0.1, fc='k', ec='k')

ax.annotate(r'$\beta + \nu_2$', xy=(0.1,.6), size=25)
ax.annotate(r'$\alpha + \nu_1$', xy=(0.3,.1), size=25)
ax.annotate(r'$O$', xy=(-0.3,.1), size=25)
ax.annotate('-3', xy=(-2.9,.1), size=25)
ax.annotate('-3', xy=(.1,-2.9), size=25)

x=[-3,-3]
y=[-3.5,.5]
ax.plot(x,y, '--', lw=3,color='black')

x=[-3.5,.5]
y=[-3,-3]
ax.plot(x,y, '--', lw=3,color='black')

x=[-3.375,0.375]
y=[0.375,-3.375]
ax.plot(x,y, '--', lw=3,color='black')

t=np.arange(-3,0.1,0.1)
upper_bound=-3-t
ax.fill_between(t, -3, upper_bound, facecolor='grey', alpha='0.5')

plt.show()
fig.savefig('bias_choice.pdf')