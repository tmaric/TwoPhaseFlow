# notebookPlotting.py

# This module contains functions used to plot data in jupyter notebooks.

import matplotlib.pyplot as plt
from matplotlib import lines
from matplotlib.lines import Line2D
from matplotlib import rcParams
import numpy as np
import glob

import dataAgglomeration as datglom

def plot_dframe(dFrame, dFrameAgglomerator, title="", plotDict = {}, ncol=2):
    """Plot the agglomerated multidimensional DataFrame using mathematical symbols mapped
    to column names with plotDict."""

    # Set line and marker styles  
    lStyles = list(lines.lineStyles.keys()) 
    lStyles = [lStyle for lStyle in lStyles if lStyle not in ('None','',' ',None)]
    mStyles = list(Line2D.markers.keys())
    mStyles = [mStyle for mStyle in mStyles if mStyle not in ('None','',' ',None)] 

    # X and Y-axis column names and diagram symbols.
    xColName = plotDict["x"]
    xColSymb = plotDict["xsymb"]
    yColName = plotDict["y"]
    yColSymb = plotDict["ysymb"]
    
    indexLevels = [i for i,name in enumerate(dFrame.index.names) if 'step' not in name]
    collector = dFrameAgglomerator.data_collector
    variations = collector.valid_variations
    fig, ax = plt.subplots()
    ax.set_xlabel(xColSymb)
    ax.set_ylabel(yColSymb)

    
    ax.set_title(title)
    variationI = 0 
    logPlot = True
    for paramLine, subDframe in dFrame.groupby(level=indexLevels):
        xCol = subDframe[xColName] 
        yCol = subDframe[yColName]
        
        # Build the parameter string paramName=paramValue for the plot legend. 
        paramString = ""
        if (len(indexLevels) > 1):
            for levelI,paramValue in enumerate(paramLine):
                indexName = dFrame.index.names[levelI]
                try:
                    indexName = plotDict[indexName]
                except:
                    pass

                paramString = paramString + indexName \
                              + "=%.1e " % paramValue
        else:
            indexName = dFrame.index.names[0]
            try:
                indexName = plotDict[indexName]
            except:
                pass
            paramString = indexName + "=%d " % paramLine 


        if (np.max(np.abs(yCol)) < 1e-15):
            logPlot = False
            
        ax.plot(xCol, yCol, label="%04d " % variationI + paramString, 
                marker=mStyles[variationI % len(mStyles)], 
                linestyle=lStyles[variationI % len(lStyles)])

        variationI = variationI + 1

    if (logPlot):
        ax.set_yscale('log') 

    ax.legend(loc='upper center', bbox_to_anchor=(0.5, -0.15),fancybox=True, shadow=True, ncol=ncol) 

    plt.show()
