### Global plotting settings
import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.ticker as tck
import matplotlib.patches as patches
import matplotlib.dates as mdates


plt.rcParams["font.family"] = 'Times New Roman'
plt.rcParams["figure.figsize"] = (18.35/2.54, 18.35/2.54 * 0.618) # Default is approx width of elsarticle.cls template with option 3p (two column width) and assumed 10pt fontsize and golden ratio for height
plt.rcParams["figure.facecolor"] = 'white'
plt.rcParams["figure.dpi"] = 300
plt.rcParams["savefig.dpi"] = 300
plt.rcParams["axes.labelsize"] = 10
plt.rcParams["axes.facecolor"] = 'white'
plt.rcParams["xtick.labelsize"] = 10
plt.rcParams["ytick.labelsize"] = 10

plt.rcParams["grid.alpha"] = 0.5
plt.rcParams["grid.linestyle"] = '-'
plt.rcParams["grid.linewidth"] = 1.25

plt.rcParams["legend.fontsize"] = 'small'
plt.rcParams["legend.framealpha"] = 0.9
plt.rcParams["legend.loc"] = 'best'