import matplotlib.pyplot as plt
import numpy as np

# -------------------------
# Data
# -------------------------
methods    = ['Our Method', 'LTR Finder', 'LTR Harvest', 'LTR Detector', 'RepeatModeler']
precision  = [0.7518, 0.7949, 0.8118, 0.8025, 0.8314]
recall     = [0.9841, 0.9757, 0.9695, 0.9344, 0.9428]
f1_scores  = [0.8524, 0.8761, 0.8837, 0.8634, 0.8836]

tp = [1294682658, 1283615787, 1275474159, 1229206879, 1240295218]
fp = [427518589,  331286097,  295706172,  302549828,  251598753]
fn = [20874309,   31941180,   40082808,   86350088,   75261749] 


method_colors = [
    'mediumspringgreen',  # Our Method
    'lightskyblue',       
    'deepskyblue',        
    'dodgerblue',         
    'cornflowerblue'
]

# Colors for categories in stacked plot (TP/FN/FP)
tp_color = 'lightskyblue'
fn_color = 'deepskyblue'
fp_color = 'dodgerblue'


n = len(methods)
assert all(len(lst) == n for lst in [precision, recall, f1_scores, tp, fp, fn, method_colors]), \
    "All lists (including method_colors) must match the number of methods."

# -------------------------
# Grouped bar chart (metrics grouped)
# -------------------------
metrics = ['Precision', 'Recall', 'F1 Score']
x_metrics = np.arange(len(metrics))
width = 0.8 / n  

fig, ax = plt.subplots(figsize=(6, 5))

for i, method in enumerate(methods):
    offset = (i - (n - 1) / 2) * width
    scores = [precision[i], recall[i], f1_scores[i]]
    ax.bar(
        x_metrics + offset,
        scores,
        width,
        label=method,
        color=method_colors[i]
    )

ax.set_xticks(x_metrics)
ax.set_xticklabels(metrics)
ax.set_ylabel('Score')
ax.set_ylim(0, 1)
ax.set_title('Maize (Zea mays)')
ax.legend(bbox_to_anchor=(1.02, 1), loc='best', borderaxespad=0.)
fig.tight_layout()
fig.savefig('maize_grouped.png', dpi=300)

# -------------------------
# Stacked bar chart (TP / FN / FP)
# -------------------------
x_methods = np.arange(n)
tp_arr = np.array(tp)
fn_arr = np.array(fn)
fp_arr = np.array(fp)

fig, ax = plt.subplots(figsize=(6, 5))

ax.bar(x_methods, tp_arr, label='TP', color=tp_color)
ax.bar(x_methods, fn_arr, bottom=tp_arr, label='FN', color=fn_color)
ax.bar(x_methods, fp_arr, bottom=tp_arr + fn_arr, label='FP', color=fp_color)

ax.set_xticks(x_methods)
ax.set_xticklabels(methods, rotation=18, ha='right')
ax.set_ylabel('Number of Bases')
ax.set_title('Maize (Zea mays)')
ax.legend(bbox_to_anchor=(1.02, 1), loc='best', borderaxespad=0.)
fig.tight_layout()
fig.savefig('maize_stacked.png', dpi=300)

plt.show()
