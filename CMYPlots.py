import pandas as pd
import numpy as np
import re
from matplotlib import pyplot as plt

def plotHighDimGrids(df,
                    xAxis = 'hourOfDay',
                    yAxis = 'dayOfWeek',
                    xLabels = 'null',
                    yLabels = 'null',
                    analyses =[('count of all edges','fromDegree','count','null')],
                    channelBy = 'columns',
                    sample= 10000,
                    scale= .5,
                    dataSetName = 'null',
                    minAtZero = False,
                    cmap = 'null',
                    analysesLabel='',
                    figOut=''):
    """Loads data, maps two int columns to x & y axis, bins rows by x & y cells, and maps analyses to CMY color channels"""
    
    cleanLabel = lambda x: ' '.join(re.findall(r'[A-Z]?[a-z]+|[A-Z]+(?=[A-Z]|$)',x)).title()
    
    if type(df) is str:
        df = getData(df)
        
    df = df[np.isfinite(df[xAxis])]
    df = df[np.isfinite(df[yAxis])]
    df.loc[:,xAxis] = df.loc[:,xAxis].astype(int)
    df.loc[:,yAxis] = df.loc[:,yAxis].astype(int)
    
    xAxisOut = cleanLabel(xAxis)
    yAxisOut = cleanLabel(yAxis)

    if sample == 'null':
        sampled = df
        sampleArg = "Over %s Entries" % len(sampled)
    else:
        sampled = df.sample(n=min(sample,len(df)))
        sampleArg = "for Sample of Size %s" % sample
        
    if type(analyses) is str:
        analyses = analyses.lower().replace('-','_').replace(' ','_')
        if analyses =='count_all':
            analysesOut = [('Count of All Edges','fromDegree','count','null')]
        elif analyses =='count_all_by_degrees':
            analysesOut = [('Low','distance','count',"df[df['fromDegreeQ']==0]"),
             ('Medium','distance','count',"df[df['fromDegreeQ']==1]"),
             ('High','distance','count',"df[df['fromDegreeQ']==2]")]
        else:
            analysisArgs = analyses.split('_')
            analysisArgs = [i for i in analysisArgs if i.lower() not in {'of','by'}]
            metric = analysisArgs[0]
            attribute = df.columns[[i.lower() for i in df.columns].index(analysisArgs[1])]
            if len(analysisArgs) == 3:
                separator = df.columns[[i.lower() for i in df.columns].index(analysisArgs[2])]
                labels = [['Low',0],['Medium',1],['High',2]]
                analysesOut = [(i,attribute,metric,"df[df['%s']==%s]" % (separator,j)) for [i,j] in labels]
                analysesLabel = cleanLabel(separator)
            elif len(analysisArgs) == 2:
                analysesOut = [('%s of all %s edges' % (cleanLabel(attribute),channel),attribute,metric,"null")]
        for i,j in enumerate(analysesOut):
            print("Analysis %s: %s" % (i,j))
        analyses = analysesOut
            
        
    xVals = range(df[xAxis].min(),df[xAxis].max()+1)
    yVals = range(df[yAxis].min(),df[yAxis].max()+1)
    nX = len(xVals)
    nY = len(yVals)
    minX = min(xVals)
    minY = min(yVals)
    
    nAnalyses = len(analyses)
    if nAnalyses > 3:
        print("Error, too many analysis channels for current algorithm")
        analyses = analyses[:3]
        nAnalyses = 3
    if nAnalyses != 1:
        cmap = 'null'
        

    xyGridCounts = pd.DataFrame(0,columns=xVals,index=yVals)
    results = dict()
    cmykValues = dict()
    
    cmykValues = np.ones([nY,nX,3])
    colorBars = []
    buffer = 3
    shape = (nX+buffer)/(nY+buffer)
    
    fig,axes = plt.subplots(nAnalyses+1,1,
                            figsize=(len(xVals)*scale,(len(yVals)+2*nAnalyses)*scale),
                            gridspec_kw={'height_ratios': [shape]+[shape/10]*nAnalyses,
                                         'width_ratios':[1/shape]})
    

    minValues = [0]*nAnalyses
    maxValues = [0]*nAnalyses
    
    for colorChannel,analysis in enumerate(analyses):
        (label,variable,metric,fx) = analysis
        results[label] = dict()
        listedValues = []
        if fx != 'null':
            fx = fx.replace('df','sampled')
            workingDf = eval(fx)
        else:
            workingDf = sampled
        for y in yVals:
            results[label][y] = dict()
            for x in xVals:
                results[label][y][x] = []
        for index, row in workingDf.iterrows():
            var = row[variable]
            if np.isfinite(var):
                results[label][row[yAxis]][row[xAxis]].append(var)
        for y in yVals:
            for x in xVals:
                entries = results[label][y][x]
                if entries == [] or entries != entries:
                    results[label][y][x] = 'null'
                else:
                    if metric == 'mean':
                        results[label][y][x] = np.mean(results[label][y][x])
                    if metric == 'median':
                        results[label][y][x] = np.median(results[label][y][x])
                    elif metric == 'max':
                        results[label][y][x] = np.max(results[label][y][x])
                    elif metric == 'sum':
                        results[label][y][x] = np.sum(results[label][y][x])
                    elif metric == 'count':
                        results[label][y][x] = len(results[label][y][x])
                    listedValues.append(results[label][y][x])
        
        if minAtZero:
            minValues[colorChannel] = 0
        else:
            minValues[colorChannel] = min(listedValues)
        maxValues[colorChannel] = max(listedValues)
        for y in yVals:
            for x in xVals:
                if results[label][y][x] != 'null':
                    cmykValues[y-minY][x-minX][colorChannel] = float(np.interp(results[label][y][x],
                                                                     [minValues[colorChannel],maxValues[colorChannel]],
                                                                     [1.,0.]))
                else:
                    cmykValues[y-minY][x-minX][colorChannel] = 1

    if nAnalyses == 1:
        for y in yVals:
            for x in xVals:
                cmykValues[y-minY][x-minX] = [cmykValues[y-minY][x-minX][0]]*3

    def addSpines(axI):
        for spine in ['bottom','top','right','left']:
            axI.spines[spine].set_color('.0')
    
    if cmap == 'null':
        axes[0].imshow(cmykValues,origin='lower')
    else:
        cmykNew = np.ones([nY,nX])
        for y in yVals:
            for x in xVals:
                cmykNew[y-minY][x-minX] = cmykValues[y-minY][x-minX][0]
        cmykValues = cmykNew
                
        axes[0].imshow(cmykValues,origin='lower',cmap=cmap)
        
    axes[0].set_xlabel(xAxisOut)
    axes[0].set_ylabel(yAxisOut)
    axes[0].grid(False)
    plt.rcParams["axes.edgecolor"] = "black"
    plt.rcParams["axes.linewidth"] = 1
    addSpines(axes[0])
    
    np.linspace(1, 0, 256)
    gradient = np.linspace(1, 0, 256)
    
    for analysisI,ax in zip(range(nAnalyses),axes[1:]):
        analysisName = analyses[analysisI][0]
        colors = []
        for i in range(256):
            if nAnalyses != 1:
                color = [1.,1.,1.]
                color[analysisI] = gradient[i]
            else:
                color = gradient[i]
            colors.append(color)
        span = maxValues[analysisI]-minValues[analysisI]
        if cmap == 'null':
            ax.imshow([colors],extent=[minValues[analysisI],maxValues[analysisI],0,span/15],cmap='gray')
        else:
            ax.imshow([colors],extent=[minValues[analysisI],maxValues[analysisI],0,span/15],cmap=cmap)
        ax.grid(False)
        ax.set_title(analysisName,fontsize=10)
        ax.get_yaxis().set_visible(False)

    for ax in axes:
        ax.axis("tight")
    
    #fig.subplots_adjust(hspace=3, wspace=0)
    
    def tryLabel(labels,key):
        try:
            return labels[key]
        except:
            return key
        
    plt.tight_layout()

    if xLabels != 'null':
        axes[0].set_xticklabels([tryLabel(xLabels,item.get_text()) for item in axes[0].get_xticklabels()])
    if yLabels != 'null':
        axes[0].set_yticklabels([tryLabel(yLabels,item.get_text()) for item in axes[0].get_yticklabels()])
    
    return fig, axes
