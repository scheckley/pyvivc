# Pyvivc

## Version 1.3
Stephen Checkley, October 2021.

---

## Description
A Python 3 port of the Rivivc R package for IVIVC linear level A by Aleksander Mendyk and Sebastian Polak. The package contains a numerical deconvolution method working for inequal and incompatible timepoints between impulse and response curves. A numerical convolution method is also included.

---

## Installation

Clone the repository and install with pip:

```
pip install .
```

or install with pip from the PyPi repository:

```
pip install pyvivc
```

## Pyvivc example

    from pyvivc import *
    from scipy.interpolate import pchip_interpolate as pchip
    from scipy.optimize import minimize
    import pandas as pd
    import numpy as np
    import matplotlib.pyplot as plt
    import matplotlib
    import seaborn as sns
    sns.set()
    matplotlib.use('TkAgg')


    impulse = pd.read_csv('data/impulse.csv')
    response = pd.read_csv('data/resp.csv')
    inp = pd.read_csv('data/input.csv')

    out = pyivivc(inp,impulse,response,explicit_interpolation=10,implicit_interpolation=5)

    rsquare_text = str('R squared = ') + str(round(out[0].rvalue,2))

    plt.subplot(1, 2, 1)
    plt.plot(inp['C'][0:len(inp)-1], out[1]['par'],'o', label='data')
    plt.plot(out[1]['par'], out[0].intercept + out[0].slope*out[1]['par'], 'r')
    plt.annotate(rsquare_text, (0,0.8), horizontalalignment='left', verticalalignment='top', fontsize=8) # these are the coordinates to position the label
    plt.xlabel('input data (#)')
    plt.ylabel('deconvolved input (#)')
    plt.legend()

    plt.subplot(1, 2, 2)
    plt.plot(out[1]['time'], out[1]['par'], 'r', label='deconvolution')
    plt.plot(inp['time'][0:len(inp)-1], inp['C'][0:len(inp)-1], 'o', label='data')
    plt.xlabel('Time')
    plt.ylabel('discovered input (%)')
    plt.legend()

    plt.show()


![The output should look like this](https://github.com/scheckley/pyvivc/blob/master/example_plot.png?raw=true)
---
