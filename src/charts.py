import numpy as np
import matplotlib.pyplot as plt
import math
import pylab as plt
import pandas as pd
import matplotlib.patches as patches

def MGM_chart(match,mismatch,gab):

  height = [match, mismatch, gab]
  bars = ('Match', 'Mismatch', 'Gabs')
  y_pos = np.arange(len(bars))
  plt.bar(y_pos, height,width=0.5,color=['purple', 'cyan','yellow'])
  plt.xticks(y_pos, bars)
  plt.show()
  

def plot_fastq_qualities(qualities, ax=None):
    '''qualities is a list represent the 4th line in fastq  '''
    scores=[]
    score=[]
    for i in range (0,len(qualities)):
      for j in range(0,len(qualities[i])):
        score.append(ord(qualities[i][j])-33)
      scores.append(score)
      score=[]
    #print(scores)

    df = pd.DataFrame(scores)
    l = len(df.T)+1

    if ax==None:
        f,ax=plt.subplots(figsize=(20,15))
    rect = patches.Rectangle((0,0),l,20,linewidth=0,facecolor='r',alpha=.4)
    ax.add_patch(rect)
    rect = patches.Rectangle((0,20),l,8,linewidth=0,facecolor='yellow',alpha=.4)
    ax.add_patch(rect)
    rect = patches.Rectangle((0,28),l,12,linewidth=0,facecolor='g',alpha=.4)
    ax.add_patch(rect)
    df.mean().plot(ax=ax,c='black')
    boxprops = dict(linestyle='-', linewidth=1, color='black')
    df.plot(kind='box', ax=ax, grid=False, showfliers=False,
            color=dict(boxes='black',whiskers='black')  )
    ax.set_xticks(np.arange(0, l, 5))
    ax.set_xticklabels(np.arange(0, l, 5))
    ax.set_xlabel('position(bp)')
    ax.set_xlim((0,l))
    ax.set_ylim((0,40))
    ax.set_title('per base sequence quality')
    plt.show()  
    
    
if __name__ == '__main__':
    
    #MGM_chart(10,5,10)
    plot_fastq_qualities(['CFFFEFHHH', 'ACGHIJIHA'])