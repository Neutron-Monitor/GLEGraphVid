#!/usr/bin/python3

"""
===========================================================================
# Make_GLEAlarm
# Script to analyze Neutron Monitor rates for Ground Level Enhacements (GLE)
# events and email alerts
#
# Auhors:
# Brian Lucas
# Pierre-Simon Mangeard
#
# Versions:
# 1.0.0 Initial version adapted from Make_GLEAlarm
"""
import glob
from datetime import datetime, timedelta, timezone, date, time
import time
import numpy as np
import pandas as pd

import math
import csv
import json
import sys
import getopt

import os.path
from os import path

import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.dates as mdates
from matplotlib.ticker import ScalarFormatter
from matplotlib.cbook import get_sample_data
import matplotlib.gridspec as gridspec

import smtplib
from email.message import EmailMessage

#######################################
###Set Latex font for figures
#######################################
mpl.rcParams.update({
	'font.size':17,
	'text.usetex': False,
	'font.family': 'stixgeneral',
	'mathtext.fontset': 'stix',
})



__author__      = 'Brian Lucas'
__credits__ = ['Brian Lucas','Pierre-Simon Mangeard']
__email__ = 'lucasb@udel.edu'

pd.options.mode.chained_assignment = None  # default='warn'


def main(argv):
   start_exetime = time.time()

   ########################
   ### DEFINE VARIABLES
   ########################
   Inpath = '.'      #input path
   Outpath = '.'     #output path
   urlalarm="./GLE_Alarm.png" #DEBUG
   initMinutes=30
   frameNum = 0

   ########################
   ### ARGUMENTS
   ########################

   strinfo='GLEGraphVid.py: options:\n'
   strinfo=strinfo+'-i <input path>\n'
   strinfo=strinfo+'-o <output path> (output path is the same as input path if not given)\n'
   strinfo=strinfo+'-r <replay day> (in valid ISO 8601 format like YYYY-MM-DD)\n'
   strinfo=strinfo+'-s <minutes>\n'
   strinfo=strinfo+'-e <minutes>\n'

   try:
      opts, args = getopt.getopt(argv,"hr:s:e:i:o:")
   except getopt.GetoptError:
      print(strinfo)
      sys.exit(2)

   for opt, arg in opts:
      if opt == '-h':
         print(strinfo)
         sys.exit()
      elif opt in ("-i"):
         Inpath = arg      #input path
         Outpath = arg     #output path
      elif opt in ("-o"):
         Outpath = arg     #output path
      elif opt in ("-r"):
         startDay = date.fromisoformat(arg) #set the start day
         print("{0:s} interpreted Replay Day {1:s}".format(arg, startDay.strftime('%D')))
      elif opt in ("-s"):
         startMinutes = int(arg) #set the start min in day
      elif opt in ("-e"):
         endMinutes = int(arg) #set the end min in day


   if len(opts) <  1:
      print('For information: GLEGraphVid.py -h')
      sys.exit(2)

   ########################
   #Data frame
   ########################

   ########################
   ### Stations to read
   ########################

   nm=      ['in','fs','pe','na','ne','th','sp','sp','mc','jb']
   nmdbtag= ['INVK','FSMT','PWNK','NAIN','NEWK','THUL','SOPO','SOPB','MCMU','JBGO']
   Labels=  ['Inuvik','Fort Smith','Peawanuck','Nain','Newark','Thule','South Pole','$^{\dagger}$South Pole - bare','McMurdo','Jang Bogo']
   InAlert= [1       ,1           ,1          ,1     ,0       ,1      ,1                       ,0                  ,1        ,0          ]
   sFact= ['','',' *2','',' *2','',' /2','','','']
   Fact=  [1.,1.,2.   ,1.,2.    ,1.,0.5 ,1.,1.,1.]

   #History from Makejson_ql.py
   #History= [0.597135,(0.598/0.94696),1.35333,0.59686,0.54518,0.57732*0.6,0.705,0.52308,0.52308]
   History= [0.597135,0.598,1.35333,0.59686,0.54518,0.57732*0.6,0.52308,0.52308,0.705]
   History=np.array(History)
   History=History/0.6
   #Number of stations
   N= len(nm)-1
   #if isReplay: N-=1 #DEBUG
   Notused='$^{\dagger}$Not used'

   df = pd.read_csv('{0:s}/GLE_Day_{1:s}.txt'.format(
                     Outpath,startDay.strftime("%Y%m%d")),
                     sep=',',date_format='%y/%m/%d %H:%M:%S',
                     index_col=0)
   df=df.drop(columns=['Time.1'])
   print(df)  #DEBUG
   print(df.info(verbose=True, show_counts=True))  #DEBUG

   startTime=datetime.combine(startDay, datetime.min.time())+timedelta(minutes=startMinutes)
   print(startTime)  #DEBUG
   endTime=datetime.combine(startDay, datetime.min.time())+timedelta(minutes=endMinutes)
   print('Graph x-axis from',startTime, 'to', endTime)

   df = df[startTime:endTime]
   print(df)  #DEBUG
   print(df.info(verbose=True, show_counts=True))  #DEBUG
   # sys.exit() #DEBUG




   ########################
   #Define Status
   ########################
   #0: Quiet
   #1: Watch
   #2: Warning
   #>=3: Alert

   Status=['Quiet','Watch','Warning','Alert']
   Statuscol=['gray','blue','orange','red']


      #print(df)  #DEBUG
      #print(df.info(verbose=True, show_counts=True))  #DEBUG







      # for i in range(N):
      #   df=df.drop(columns=[nmdbtag[i]+'F'])



   ########################
   ### Display
   ########################



   fontsize=17

   ymaxT=30000.
   yminT=0.
   ymaxI=10.
   yminI=-5.
   
   for i in range(N):
      if df[nmdbtag[i]+'T'].max()>ymaxT: ymaxT=df[nmdbtag[i]+'T'].max()
      if df[nmdbtag[i]+'T'].min()<yminT: yminT=df[nmdbtag[i]+'T'].min()
      if df[nmdbtag[i]+'Ith'].isnull().values.any(): pass #print('Null values in {0:s} during plotting'.format(nmdbtag[i]+'Ith'))
      else:
         if 100.*(df[nmdbtag[i]+'Ith'].max()-1.)>ymaxI: ymaxI=100.*(df[nmdbtag[i]+'Ith'].max()-1.)
         if 100.*(df[nmdbtag[i]+'Ith'].min()-1.)<yminI: yminI=100.*(df[nmdbtag[i]+'Ith'].min()-1.)


   LastStatus=0
   for r in range(initMinutes, endMinutes- startMinutes) :
         
      dfCur=df.iloc[0:r]
      fig=plt.figure(figsize=(14, 11), dpi=80)

      axesT = fig.add_subplot(5,1,(2,3))

      for i in range(N):
         # axesT.plot(df['Time'],Fact[i]*df[nmdbtag[i]+'T'],'-',linewidth=0.8,label=r'{0:s} {1:s}'.format(Labels[i],sFact[i]))
         axesT.plot(dfCur.index.values,Fact[i]*dfCur[nmdbtag[i]+'T'],'-',linewidth=0.8,label='{0:s} {1:s}'.format(Labels[i],sFact[i]))
         # if df[nmdbtag[i]+'T'].max()>ymax: ymax=df[nmdbtag[i]+'T'].max()
         # if df[nmdbtag[i]+'T'].min()<ymin: ymin=df[nmdbtag[i]+'T'].min()

      plt.tick_params(axis='x', which='major', labelsize=0,direction='in',length=6)
      plt.tick_params(axis='y', which='major', labelsize=fontsize,direction='in',length=6)
      plt.tick_params(axis='y', which='major', labelsize=fontsize,direction='in',length=6)

      axesT.xaxis.set_major_locator(mdates.HourLocator(interval=1))
      axesT.xaxis.set_major_formatter(mdates.DateFormatter('%m/%d\n%H:%M'))
      # axesT.set_xlim(now - timedelta(hours=12),now)
      axesT.set_xlim(startTime,endTime)
      axesT.set_ylim(yminT,ymaxT-1)
      axesT.set_ylabel('Rate [count / minute]\n3-min moving average',fontsize=fontsize+1)

      plt.grid(axis='both',which='major',linewidth=0.5,linestyle=':',color='gray')

      plt.legend(bbox_to_anchor=(1.01,0.525), loc="center left", borderaxespad=0,
               fontsize=fontsize,labelspacing=1,frameon=False)
      plt.title('GLE Alarm from Bartol Neutron Monitors - {0:s} UT'.format(dfCur.index.values[-1]),fontsize=fontsize)
      # plt.title('GLE Alarm from Bartol Neutron Monitors - {0:s} {1:s} {2:s} UT'.format(df.index.values[r],r'$@$',r.strftime("%H:%M:%S")),fontsize=fontsize)

      axes = fig.add_subplot(5,1,(4,5),sharex=axesT)
      

      for i in range(N):
         axes.plot(dfCur.index.values,100.*(dfCur[nmdbtag[i]+'Ith']-1.),'-',linewidth=0.8,label='{0:s}'.format(Labels[i]))


      plt.tick_params(axis='x', which='major', labelsize=fontsize+1,direction='in',length=6)
      plt.tick_params(axis='y', which='major', labelsize=fontsize,direction='in',length=6)
      # plt.tick_params(axis='y', which='major', labelsize=fontsize,direction='in',length=6)

      axes.xaxis.set_major_locator(mdates.HourLocator(interval=1))
      axes.xaxis.set_major_formatter(mdates.DateFormatter('%m/%d\n%H:%M'))
      ###axes.set_xlim(now - timedelta(hours=8),now)
      axes.set_ylim(yminI,ymaxI-0.001)
      axes.set_ylabel('Rate increase [%]\n(3-min tr. moving average)',fontsize=fontsize+1)
      axes.axhline(y=dfCur.iloc[-1]['Status'],linewidth=0.5,linestyle='-',color='red')
      # axes.axhline(y=Level,linewidth=0.5,linestyle='-',color='red')

      plt.grid(axis='both',which='major',linewidth=0.5,linestyle=':',color='gray')
      plt.legend(bbox_to_anchor=(1.01,0.525), loc="center left", borderaxespad=0,
               fontsize=fontsize,labelspacing=1,frameon=False)
      #plt.title('GLE Alarm from Bartol Neutron Monitors - Last update: {0:s} {1:s} {2:s} UT'.format(now.strftime("%Y-%m-%d"),r'$@$',now.strftime("%H:%M:%S")),fontsize=fontsize)
      plt.title('GLE Alarm from Bartol Neutron Monitors - Last update: {0:s} UT'.format(dfCur.index.values[-1]),fontsize=fontsize)

      axes.text(axes.get_xlim()[1] + 0.19*(axes.get_xlim()[1] -axes.get_xlim()[0] ) ,
               axes.get_ylim()[0]+ 0.0*(axes.get_ylim()[1] -axes.get_ylim()[0] ),
               Notused, horizontalalignment='left', fontsize=fontsize-2,zorder=10)

      
      axesal = fig.add_subplot(5,1,(1,1),sharex=axes)
      dfStatus = dfCur[df['Status']==3]
      # print(dfStatus) #DEBUG
      axesal.plot(dfStatus.index.values,3*np.ones(len(dfStatus)),'o',color='red',label=Status[3])
      dfStatus = dfCur[df['Status']==2]
      axesal.plot(dfStatus.index.values,2*np.ones(len(dfStatus)),'o',color='orange',label=Status[2])
      dfStatus = dfCur[df['Status']==1]
      axesal.plot(dfStatus.index.values,1*np.ones(len(dfStatus)),'o',color='blue',label=Status[1])
      dfStatus = dfCur[df['Status']==0]
      axesal.plot(dfStatus.index.values,0*np.ones(len(dfStatus)),'o',color='gray',label=Status[0])
      axesal.set_ylabel('Alarm Level',fontsize=fontsize+1)
      # axesal.plot(df[df['Status']==3]['Time'],3*np.ones(len(df[df['Status']==3])),'o',color='red',label=Status[3])
      # axesal.plot(df[df['Status']==2]['Time'],2*np.ones(len(df[df['Status']==2])),'o',color='orange',label=Status[2])
      # axesal.plot(df[df['Status']==1]['Time'],1.*np.ones(len(df[df['Status']==1])),'o',color='blue',label=Status[1])
      # axesal.plot(df[df['Status']==0]['Time'],0*np.ones(len(df[df['Status']==0])),'o',color='gray',label=Status[0])
      # axesal.set_ylabel('Alarm Level',fontsize=fontsize+1)

      plt.tick_params(axis='x', which='major', labelsize=0,direction='in',length=6)
      plt.tick_params(axis='x', which='minor', labelsize=0,direction='in',length=3)
      plt.tick_params(axis='y', which='major', labelsize=fontsize,direction='in',length=6)
      plt.yticks(np.arange(0, 4, 1.0))

      axesal.xaxis.set_major_locator(mdates.HourLocator(interval=1))
      axesal.xaxis.set_major_formatter(mdates.DateFormatter('%m/%d\n%H:%M'))
      axesal.set_ylim(0,3.75)
      plt.grid(axis='x',which='major',linewidth=0.5,linestyle=':',color='gray')
      plt.legend(bbox_to_anchor=(1.01,0.48), loc="center left", borderaxespad=0,
               fontsize=fontsize,labelspacing=1,frameon=False)

      axesal.text(axesal.get_xlim()[0] + 0.02*(axesal.get_xlim()[1] -axesal.get_xlim()[0] ) ,
               axesal.get_ylim()[1]- 0.12*(axesal.get_ylim()[1] -axesal.get_ylim()[0] ),
               "Current: ", horizontalalignment='left', fontsize=fontsize+1,zorder=10)

      axesal.text(axesal.get_xlim()[0] + 0.11*(axesal.get_xlim()[1] -axesal.get_xlim()[0] ) ,
               axesal.get_ylim()[1]- 0.12*(axesal.get_ylim()[1] -axesal.get_ylim()[0] ),
               "{0:s}".format(Status[int(LastStatus)]), horizontalalignment='left',color=Statuscol[int(LastStatus)], fontsize=fontsize+1,zorder=10)


      axesal.text(axesal.get_xlim()[0], axesal.get_ylim()[1]+ 0.017*(axesal.get_ylim()[1] -axesal.get_ylim()[0] ),
            'Last update: {0:s} UT'.format(dfCur.index.values[-1]),fontsize=fontsize+1,horizontalalignment='left')


      plt.subplots_adjust(left=0.1, bottom=0.06, right=0.8, top=0.95, wspace=0, hspace=0.00)
   
      # fig.savefig('{0:s}/GLE_Alarm.png'.format(Outpath))
      fig.savefig('{0:s}/{1:s}/{2:04d}.png'.format(Outpath, startTime.strftime("%Y%m%d"), frameNum))
      frameNum+=1
      LastStatus=int(dfCur.iloc[-1]['Status'])
      plt.clf()


   #for i in range(N-1):
   #   df=df.drop(columns=[nmdbtag[i+1]+'T'])
   #   df=df.drop(columns=[nmdbtag[i+1]+'Ith'])
   #   df=df.drop(columns=[nmdbtag[i+1]])

   #print(df.iloc[-10:])
   #print("--- %s seconds ---" % (time.time() - start_exetime))

   # plt.show()



   # print(df.info(verbose=True, show_counts=True))  #DEBUG
   #print(df[-2:])  #DEBUG
   # print(archive_data.info(verbose=True, show_counts=True))  #DEBUG

   # sys.exit() #DEBUG


if __name__ == "__main__":
   main(sys.argv[1:])



#END
