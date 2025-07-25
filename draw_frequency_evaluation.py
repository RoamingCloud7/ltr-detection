import matplotlib.pyplot as plt
import numpy as np

# ---------------------
# Data
# ---------------------
classes = ['High-copy', 'Medium-copy', 'Low-copy']
x = np.arange(len(classes))

recall_mine           = [0.989270, 0.902526, 0.444876]
recall_finder         = [0.982564, 0.864291, 0.575545]
recall_harvest        = [0.978324, 0.827887, 0.276199]
recall_detector       = [0.933618, 0.947868, 0.821713]
recall_repeatmodeler  = [0.949455, 0.836861, 0.258747]
 

# ---------------------
# Plot
# ---------------------
fig, ax = plt.subplots(figsize=(6, 4.5))

ax.plot(x, recall_mine,           marker='o', label='Our Method',       color='mediumspringgreen')
ax.plot(x, recall_finder,        marker='o', label='LTR Finder',         color='lightskyblue') 
ax.plot(x, recall_harvest,        marker='o', label='LTR Harvest',      color='deepskyblue')
ax.plot(x, recall_detector,       marker='o', label='LTR Detector',     color='dodgerblue')
ax.plot(x, recall_repeatmodeler,  marker='o', label='RepeatModeler',    color='cornflowerblue')


# ---------------------
# Styling
# ---------------------
ax.set_xticks(x)
ax.set_xticklabels(classes)
ax.set_ylim(0, 1)
ax.set_ylabel('Recall')
ax.set_title('Maize (Zea mays)')
ax.grid(True, axis='y')


ax.legend(loc='best', frameon=True, framealpha=0.85)

fig.tight_layout()
fig.savefig('maize_line.png', dpi=300)
plt.show()
