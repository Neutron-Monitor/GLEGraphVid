#!/usr/bin/python3

"""
===========================================================================
# GLEGraphVid.py
# Script to make Graphs of GLE to compile into a video
#
# Auhors:
# Brian Lucas
# Pierre-Simon Mangeard
#
# Versions:
# 1.0.0 Initial version adapted from Make_GLEAlarm
# 1.1.0 Add GOES data
# 1.2.0 Add vertical alarm lines
# 1.3.0 Add vertical baselines
# 1.4.0 Add margins to limits
# 1.5.0 Change margins using with 0 ymin as special case
# 1.6.0 Change in stations
# 1.7.0 Change in stations, handling of not in alert and tickmarks
# 1.8.0 Added different NM network awareness
# 1.9.0 Step plot for network aware. Option to exclude rates plot
# 1.10.0 Interpret new GOES data format
# 1.11.0 Formatting changes for GLE77
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
# from matplotlib.ticker import ScalarFormatter
import matplotlib.ticker as mticker
from matplotlib.cbook import get_sample_data
import matplotlib.gridspec as gridspec

import smtplib
from email.message import EmailMessage
from cycler import cycler



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
   # urlalarm='./GLE_Alarm.png' #DEBUG
   # initMinutes=90
   initMinutes=1
   frameNum = 0
   fileGOESProton =''
   fileGOESXray ='' 
   showBaselines = False
   alarmLineGPShow = True
   xTickMajorHours = 1

   ########################
   ### ARGUMENTS
   ########################

   strinfo='GLEGraphVid.py: options:\n'
   strinfo=strinfo+'-i <input path>\n'
   strinfo=strinfo+'-o <output path> (output path is the same as input path if not given)\n'
   strinfo=strinfo+'-r <replay day> (in valid ISO 8601 format like YYYY-MM-DD)\n'
   strinfo=strinfo+'-s <minutes>\n'
   strinfo=strinfo+'-e <minutes>\n'
   strinfo=strinfo+'-p <GOES Proton file>\n'
   strinfo=strinfo+'-x <GOES X-ray file>\n'
   strinfo=strinfo+'-b (show baseline)\n'

   try:
      opts, args = getopt.getopt(argv,"hr:s:e:i:o:p:x:b")
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
      elif opt in ("-p"):
         fileGOESProton = arg     #GOES Proton file
      elif opt in ("-x"):
         fileGOESXray = arg     #GOES X-ray file
      elif opt in ("-b"):
         showBaselines = True     #graph baselines


   if len(opts) <  1:
      print('For information: GLEGraphVid.py -h')
      sys.exit(2)

   ########################
   #Data frame
   ########################

   ########################
   ### Stations to read
   ########################

   nm=      ['in','fs','pe','na','ne','th','sp','sp','mc','jb','','','','','','']
   nmdbtag= ['INVK','FSMT','PWNK','NAIN','NEWK','THUL','SOPO','SOPB','MCMU','JBGO','MWSN','CVAN','DRHM','HLE1','LDVL','MTWS']
   Labels=  ['Inuvik','Fort Smith','Peawanuck','Nain','Newark','Thule','South Pole','South Pole - bare','McMurdo','Jang Bogo','Mawson','ChangVan','Durham','Haleakalā','Leadville','Mt. Washington']
   InAlert= [1       ,1           ,1          ,1     ,1       ,1      ,1           ,0                  ,0        ,0          ,0       ,0,0,0,0,0]
   sFact= ['','',' *2','',' *2','',' /2','','','','','','','','','']
   Fact=  [1.,1.,2.   ,1.,2.    ,1.,0.5 ,1.,1.,1.,1.,1.,1.,1.,1.,1.]

   bartolFlags = ['INVKF','FSMTF','PWNKF','NAINF','NEWKF','THULF', 'SOPOF']
   extendedFlags = ['SOPBF','DRHMF','HLE1F','LDVLF','MTWSF']
   intlFlags = list(map(lambda x: x + 'F', nmdbtag))
   intlFlags = list(filter(lambda x: x not in bartolFlags, intlFlags))
   intlFlags = list(filter(lambda x: x not in extendedFlags, intlFlags))
   print(intlFlags)  #DEBUG

   #History from Makejson_ql.py
   #History= [0.597135,(0.598/0.94696),1.35333,0.59686,0.54518,0.57732*0.6,0.705,0.52308,0.52308]
   History= [0.597135,0.598,1.35333,0.59686,0.54518,0.57732*0.6,0.52308,0.52308,0.705,0.6,0.6,0.6,0.6,0.6,0.6]
   History=np.array(History)
   History=History/0.6

   lenGP=0
   lenGX=0
   limMargin = 0.01

   networkAwareAlert = True
   ratePlot = False

   # df = pd.read_csv('{0:s}/GLE_Day_{1:s}.txt'.format(
   df = pd.read_csv('{0:s}/GLE_Day_{1:s}.csv'.format(
                     Outpath,startDay.strftime("%Y%m%d")),
                     sep=',',date_format='%y/%m/%d %H:%M:%S',
                     index_col=0)
   df=df.drop(columns=['Time.1'])
   # print(df)  #DEBUG
   # print(df.info(verbose=True, show_counts=True))  #DEBUG


   #Number of stations
   N= len(nm)
   Nall = N
   Notused='$^{\dagger}$Not used'
   i=1
   while i < N:
      if not InAlert[i]:
         nm.append(nm.pop(i)    )
         nmdbtag.append(nmdbtag.pop(i))
         Labels.append('$^{\dagger}$'+Labels.pop(i) )
         InAlert.append(InAlert.pop(i))
         sFact.append(sFact.pop(i))
         Fact.append(Fact.pop(i))
         N-=1
      else : i+=1
   colorsPlot = plt.colormaps["tab20"](np.linspace(0, 1, Nall))
   plt.rcParams["axes.prop_cycle"] = cycler(color=colorsPlot)
   # print(nm)
   # print(nmdbtag)
   # print(Labels)
   # print(InAlert)
   # print(sFact)
   # print(Fact)
   # sys.exit() #DEBUG


   startTime=datetime.combine(startDay, datetime.min.time())+timedelta(minutes=startMinutes)
   # print(startTime)  #DEBUG
   endTime=datetime.combine(startDay, datetime.min.time())+timedelta(minutes=endMinutes)
   print('Graph x-axis from',startTime, 'to', endTime)

   df = df[startTime:endTime]
   # df = df[startTime:(endTime+timedelta(minutes=1))]
   print(df)  #DEBUG
   # print(df.info(verbose=True, show_counts=True))  #DEBUG

   alarmLines = [df[df['Status']>=1].index[0],df[df['Status']>=2].index[0],df[df['Status']>=3].index[0]]
   # alarmLines = [df[df['Status']>=1].index[0],df[df['Status']>=2].index[0],df[df['Status']>=3].index[0]]
   
   # print(alarmLines)  #DEBUG
   alarmColors = ['blue','orange','red']
   # sys.exit() #DEBUG

   lineGP=453
   if os.path.isfile(fileGOESProton):
      if ('GOES' in fileGOESProton):
         dfGPin = pd.read_csv(fileGOESProton,
                        sep=',',date_format='%y-%m-%d %H:%M:%S.%f',parse_dates=['time_tag'],
                        index_col=0, na_values=np.nan)
                        # usecols=[26,30])
                        # usecols=['p3_flux_ic','p7_flux_ic'])
         dfGPin.index = pd.to_datetime(dfGPin.index).tz_localize(None)
         # print(dfGPin.info(verbose=True, show_counts=True))  #DEBUG
         # dfGP= pd.DataFrame(index=dfGPin.loc[~dfGPin.index.duplicated(keep='first'), :].index.values)
         dfGP = dfGPin['>=10 MeV'==dfGPin['energy']]
         dfGP = dfGP.drop(['satellite','energy'], axis=1)
         dfGP=dfGP.loc[~dfGP.index.duplicated(keep='first'), :]
         dfGP.columns=['p3_flux_ic']

         dfGPp7 = dfGPin['>=100 MeV'==dfGPin['energy']]
         dfGPp7=dfGPp7.loc[~dfGPp7.index.duplicated(keep='first'), :]

         dfGP['p7_flux_ic'] = dfGPp7['flux']

         print(dfGP)  #DEBUG
         

         # print(dfGP.info(verbose=True, show_counts=True))  #DEBUG
         # sys.exit() #DEBUG
      
      else:
         lineGPoff = -5
         with open(fileGOESProton,'r') as offFile:
            for l in offFile.readlines()[lineGP+lineGPoff-1:lineGP+5]:
               if 'data' in l:
                  lineGP+=lineGPoff
                  # print(lineGP) #DEBUG
                  break
               lineGPoff+=1

         dfGP = pd.read_csv(fileGOESProton,
                        sep=',',date_format='%y-%m-%d %H:%M:%S.%f',parse_dates=['time_tag'],
                        index_col=0, skiprows=lineGP, na_values=np.nan)
                        # usecols=[26,30])
                        # usecols=['p3_flux_ic','p7_flux_ic'])
         dfGP.index = pd.to_datetime(dfGP.index)

         for c in ['e1_flux_ic','e2_flux_ic','e3_flux_ic','p1_flux','p2_flux','p3_flux','p4_flux','p5_flux','p6_flux','p7_flux','a1_flux','a2_flux','a3_flux','a4_flux','a5_flux','a6_flux','p1_flux_c','p2_flux_c','p3_flux_c','p4_flux_c','p5_flux_c','p6_flux_c','p7_flux_c','p1_flux_ic','p2_flux_ic','p4_flux_ic','p5_flux_ic','p6_flux_ic']:
            dfGP=dfGP.drop(columns=c) #drop all but p3_flux_ic, p7_flux_ic
         dfGP = dfGP[startTime:endTime]
      lenGP=len(dfGP)
      # print(dfGP[dfGP['p3_flux_ic']>=10].index[0]) #DEBUG
      # print(dfGP[dfGP['p3_flux_ic']>=100].index[0]) #DEBUG
      # print(dfGP[dfGP['p7_flux_ic']>=1].index[0]) #DEBUG
      # print(dfGP[dfGP['p7_flux_ic']>=10].index[0]) #DEBUG
      # print(alarmLines) #DEBUG
      # print(dfGP)  #DEBUG
      # print(dfGP.info(verbose=True, show_counts=True))  #DEBUG
      # sys.exit() #DEBUG

   lineX=117
   if os.path.isfile(fileGOESXray):
      if ('GOES' in fileGOESXray):
         dfGXin = pd.read_csv(fileGOESXray,
                        sep=',',date_format='%y-%m-%d %H:%M:%S.%f',parse_dates=['time_tag'],
                        index_col=0, na_values=np.nan)
                        # usecols=[26,30])
                        # usecols=['p3_flux_ic','p7_flux_ic'])
         dfGXin.index = pd.to_datetime(dfGXin.index).tz_localize(None)
         # print(dfGPin.info(verbose=True, show_counts=True))  #DEBUG
         # dfGP= pd.DataFrame(index=dfGPin.loc[~dfGPin.index.duplicated(keep='first'), :].index.values)
         dfGX = dfGXin['0.05-0.4nm'==dfGXin['energy']]
         dfGX=dfGX.loc[~dfGX.index.duplicated(keep='first'), :]
         dfGX = dfGX.drop(['satellite','observed_flux','electron_correction','electron_contaminaton','energy'], axis=1)
         dfGX.columns=['xs']
         # print(dfGX)  #DEBUG

         dfGXl = dfGXin['0.1-0.8nm'==dfGXin['energy']]
         dfGXl=dfGXl.loc[~dfGXl.index.duplicated(keep='first'), :]
         # print(dfGXl)  #DEBUG

         dfGX['xl'] = dfGXl['flux']
         dfGX=dfGX.mask(0.0==dfGX)
         # print(dfGX)  #DEBUG

         

         # print(dfGP.info(verbose=True, show_counts=True))  #DEBUG
         # sys.exit() #DEBUG
      else:
         lineXoff = -5
         with open(fileGOESXray,'r') as offFile:
            for l in offFile.readlines()[lineX+lineXoff-1:lineX+5]:
               if 'data' in l:
                  lineX+=lineXoff
                  # print(lineX) #DEBUG
                  break
               lineXoff+=1
         # sys.exit() #DEBUG

         dfGX = pd.read_csv(fileGOESXray,
                        sep=',',date_format='%y-%m-%d %H:%M:%S.%f',parse_dates=['time_tag'],
                        index_col=0, skiprows=lineX, na_values=np.nan)
                        # usecols=[26,30])
                        # usecols=['p3_flux_ic','p7_flux_ic'])
         dfGX.index = pd.to_datetime(dfGX.index)

         # for c in ['e1_flux_ic','e2_flux_ic','e3_flux_ic','p1_flux','p2_flux','p3_flux','p4_flux','p5_flux','p6_flux','p7_flux','a1_flux','a2_flux','a3_flux','a4_flux','a5_flux','a6_flux','p1_flux_c','p2_flux_c','p3_flux_c','p4_flux_c','p5_flux_c','p6_flux_c','p7_flux_c','p1_flux_ic','p2_flux_ic','p4_flux_ic','p5_flux_ic','p6_flux_ic']:
         #    dfGP=dfGP.drop(columns=c) #drop all but p3_flux_ic, p7_flux_ic
         dfGX = dfGX[startTime:endTime]
      lenGX=len(dfGX)
      # print(dfGX)  #DEBUG
      # print(dfGX.info(verbose=True, show_counts=True))  #DEBUG
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

   pT=0
   if (ratePlot):
      ymaxT=5000.
      yminT=0.
      pT=2
   ymaxI=170.
   yminI=0.

   # for i in range(N):
   for i in range(Nall):
      if (ratePlot):
         if (Fact[i]*df[nmdbtag[i]+'T'].max())>ymaxT: ymaxT=(Fact[i]*df[nmdbtag[i]+'T'].max())
         if (Fact[i]*df[nmdbtag[i]+'T'].min())<yminT: yminT=(Fact[i]*df[nmdbtag[i]+'T'].min())
         # if (1.+limMargin)*(Fact[i]*df[nmdbtag[i]+'T'].max())>ymaxT: ymaxT=(1.+limMargin)*(Fact[i]*df[nmdbtag[i]+'T'].max())
         # if (1.-limMargin)*(Fact[i]*df[nmdbtag[i]+'T'].min())<yminT: yminT=(1.-limMargin)*(Fact[i]*df[nmdbtag[i]+'T'].min())
         # if ymaxT > 20000: print(nmdbtag[i])
         # if df[nmdbtag[i]+'Ith'].isnull().values.any(): pass #print('Null values in {0:s} during plotting'.format(nmdbtag[i]+'Ith'))
      if False: pass #print('Null values in {0:s} during plotting'.format(nmdbtag[i]+'Ith'))
      else:
         if 100.*(df[nmdbtag[i]+'Ith'].max()-1)>ymaxI: ymaxI=100.*(df[nmdbtag[i]+'Ith'].max()-1)
         if 100.*(df[nmdbtag[i]+'Ith'].min()-1)<yminI: yminI=100.*(df[nmdbtag[i]+'Ith'].min()-1)
      if (ratePlot):
         ydeltaT=ymaxT-yminT
         ymaxT+=(limMargin*ydeltaT)
         yminT-=(limMargin*ydeltaT)
         if 0.>yminT :
            yminT=0.
         
      ydeltaI=ymaxI-yminI
      ydeltaI/=100
      ymaxI+=(limMargin*ydeltaI)
      if 0.!=yminI :
         yminI-=(limMargin*ydeltaI)
   pG=0
   if lenGP>0:
      ymaxGP=1.
      yminGP=0.2
      pG+=1
      # print(dfGP['p3_flux_ic'].max())
      if dfGP['p3_flux_ic'].max()>ymaxGP: ymaxGP=dfGP['p3_flux_ic'].max()
      if dfGP['p3_flux_ic'].min()<yminGP: yminGP=dfGP['p3_flux_ic'].min()
      if dfGP['p7_flux_ic'].max()>ymaxGP: ymaxGP=dfGP['p7_flux_ic'].max()
      if dfGP['p7_flux_ic'].min()<yminGP: yminGP=dfGP['p7_flux_ic'].min()
      # if (1.+limMargin)*dfGP['p3_flux_ic'].max()>ymaxGP: ymaxGP=(1.+limMargin)*dfGP['p3_flux_ic'].max()
      # if (1.-limMargin)*dfGP['p3_flux_ic'].min()<yminGP: yminGP=(1.-limMargin)*dfGP['p3_flux_ic'].min()
      # if (1.+limMargin)*dfGP['p7_flux_ic'].max()>ymaxGP: ymaxGP=(1.+limMargin)*dfGP['p7_flux_ic'].max()
      # if (1.-limMargin)*dfGP['p7_flux_ic'].min()<yminGP: yminGP=(1.-limMargin)*dfGP['p7_flux_ic'].min()
      ydeltaGP=math.log10(ymaxGP)-math.log10(yminGP)
      print(ydeltaGP) #DEBUG
      ymaxGP=10**(math.log10(ymaxGP)+limMargin*ydeltaGP)
      yminGP=10**(math.log10(yminGP)-limMargin*ydeltaGP)
      if 2e-3>yminGP :
         yminGP=2e-3

   if lenGX>0:
      ymaxGX=1e-4
      yminGX=1e-4
      # if (1.+10*limMargin)*dfGX['xs'].max()>ymaxGX: ymaxGX=(10.+10*limMargin)*dfGX['xs'].max()
      # if (1.-10*limMargin)*dfGX['xs'].min()<yminGX: yminGX=(1.-10*limMargin)*dfGX['xs'].min()
      if dfGX['xl'].max()>ymaxGX: ymaxGX=dfGX['xl'].max()
      if dfGX['xl'].min()<yminGX: yminGX=dfGX['xl'].min()
      if dfGX['xs'].max()>ymaxGX: ymaxGX=dfGX['xs'].max()
      if dfGX['xs'].min()<yminGX: yminGX=dfGX['xs'].min()
      # if (1.+limMargin)*dfGX['xl'].max()>ymaxGX: ymaxGX=(1.+limMargin)*dfGX['xl'].max()
      # if (1.-limMargin)*dfGX['xl'].min()<yminGX: yminGX=(1.-limMargin)*dfGX['xl'].min()
      ydeltaGX=math.log10(ymaxGX)-math.log10(yminGX)
      ymaxGX=10**(math.log10(ymaxGX)+limMargin*ydeltaGX)
      yminGX=10**(math.log10(yminGX)-limMargin*ydeltaGX)
      # if 2e-7>yminGX :
      #    yminGX=2e-7
      pG+=1

   if networkAwareAlert :
      df['Bartol_Above']= df[bartolFlags].sum(axis=1)
      print(df['Bartol_Above'].max())  #DEBUG
      print(df['Bartol_Above'])  #DEBUG
      df['Extended_Above']= df[extendedFlags].sum(axis=1)
      print(df['Extended_Above'].max())  #DEBUG
      df['Intl_Above']= df[intlFlags].sum(axis=1)
      print(df['Intl_Above'].max())  #DEBUG
      ymaxAl = df[['Bartol_Above','Extended_Above','Intl_Above']].sum(axis=1).max()
      ymaxAl = max(ymaxAl,4)
      print("Max Flags = {0:d}".format(ymaxAl))  #DEBUG

      df.to_csv('./GLETemp.csv')  #DEBUG

   # pAll=5+pG
   pAll=3+pG+pT
   LastStatus=0
   fig=plt.figure(figsize=(14, 11), dpi=80)
   if showBaselines: baselines = [df.index[initMinutes-85],df.index[initMinutes-10] ]
   
   # print(yminT,ymaxT,yminI,ymaxI,yminGP,ymaxGP,yminGX,ymaxGX) #DEBUG
   # print(df.index[(endMinutes)-startMinutes]) #DEBUG

   # print(df) #DEBUG
   # print (df.index[0],df.index[-1]) #DEBUG
   for r in range(initMinutes, (endMinutes+1)-startMinutes) :
   # for r in range(initMinutes, endMinutes-startMinutes) :
      # print(df.index[r]) #DEBUG
   # for r in range(endMinutes-startMinutes-1, endMinutes-startMinutes) :

      # dfCur=df.iloc[0:r]
      dfCur=df.iloc[0:(r+1)]
      # print(len(dfCur)) #DEBUG
      # print(dfCur.index[-1]) #DEBUG

      if lenGP>0:
         dfGPCur = dfGP[startTime:dfCur.index.values[-1]]
         # print(dfGPCur)  #DEBUG
         # print(dfGPCur.info(verbose=True, show_counts=True))  #DEBUG
      if lenGX>0:
         dfGXCur = dfGX[startTime:dfCur.index.values[-1]]
         # print(dfGXCur)  #DEBUG
         # print(dfGXCur.info(verbose=True, show_counts=True))  #DEBUG

      # axes = fig.add_subplot(5,1,(4,5),sharex=axesT)
      axes = fig.add_subplot(pAll,1,(pAll-1,pAll))


      for i in range(N):
         axes.plot(dfCur.index.values,100.*(dfCur[nmdbtag[i]+'Ith']-1.),'-',linewidth=0.8,label='{0:s}'.format(Labels[i]))
      for i in range(N,Nall):
         if(0 < dfCur[nmdbtag[i]+'Ith'].count()): 
            axes.plot(dfCur.index.values,100.*(dfCur[nmdbtag[i]+'Ith']-1.),'-',linewidth=0.8,label='{0:s}'.format(Labels[i]))


      axes.xaxis.set_major_locator(mdates.HourLocator(interval=xTickMajorHours))
      axes.xaxis.set_minor_locator(mdates.MinuteLocator(byminute=range(0,xTickMajorHours*60,xTickMajorHours*15)))
      axes.xaxis.set_major_formatter(mdates.DateFormatter('%Y/%m/%d\n%H:%M'))
      

      # plt.tick_params(axis='x', which='major', labelsize=fontsize+1,direction='in',length=6)
      # plt.tick_params(axis='y', which='major', labelsize=fontsize,direction='in',length=6)
      plt.tick_params(axis='x', which='major', labelsize=fontsize+1,direction='out',length=6)
      plt.tick_params(axis='x', which='minor', labelsize=0,direction='out',length=3)
      plt.tick_params(axis='y', which='major', labelsize=fontsize,direction='in',length=6)
      # plt.tick_params(axis='y', which='minor', labelsize=0,direction='in',length=3)

      ###axes.set_xlim(now - timedelta(hours=8),now)
      # axes.set_xlim(startTime,endTime)
      axes.set_xlim(df.index[0],df.index[-1])
      # if (0==yminI):
      #    axes.set_ylim(yminI,(1+limMargin)*ymaxI)
      # else:
      #    axes.set_ymargin(limMargin)
         # axes.set_ylim(yminI,ymaxI)
      axes.set_ylim(yminI,ymaxI)
      # axes.set_ylim(yminI,ymaxI-0.001)
      axes.set_ylabel('Rate increase [%]\n3-min moving average',fontsize=fontsize+1)
      # axes.axhline(y=dfCur.iloc[-1]['Status'],linewidth=0.5,linestyle='-',color='red')
      # axes.axhline(y=Level,linewidth=0.5,linestyle='-',color='red')
      

      plt.grid(axis='both',which='both',linewidth=0.5,linestyle=':',color='gray')
      plt.legend(bbox_to_anchor=(1.01,0.525), loc="center left", borderaxespad=0,
               fontsize=fontsize,labelspacing=0.0,frameon=False)
      #plt.title('GLE Alarm from Bartol Neutron Monitors - Last update: {0:s} {1:s} {2:s} UT'.format(now.strftime("%Y-%m-%d"),r'$@$',now.strftime("%H:%M:%S")),fontsize=fontsize)
      # plt.title('GLE Alarm from Bartol Neutron Monitors - Last update: {0:s} UT'.format(dfCur.index.values[-1]),fontsize=fontsize)

      axes.text(axes.get_xlim()[1] + 0.19*(axes.get_xlim()[1] -axes.get_xlim()[0] ) ,
               axes.get_ylim()[0]+ 0.0*(axes.get_ylim()[1] -axes.get_ylim()[0] ),
               Notused, horizontalalignment='left', fontsize=fontsize-2,zorder=10)

      if (ratePlot):
         # axesT = fig.add_subplot(5,1,(2,3))
         axesT = fig.add_subplot(pAll,1,(pAll-3,pAll-2),sharex=axes)
         # axesT = fig.add_subplot(pAll,1,(pAll-2,pAll-1))

         for i in range(N):
            # axesT.plot(df['Time'],Fact[i]*df[nmdbtag[i]+'T'],'-',linewidth=0.8,label=r'{0:s} {1:s}'.format(Labels[i],sFact[i]))
            axesT.plot(dfCur.index.values,Fact[i]*dfCur[nmdbtag[i]+'T'],'-',linewidth=0.8,label='{0:s} {1:s}'.format(Labels[i],sFact[i]))
            # if df[nmdbtag[i]+'T'].max()>ymax: ymax=df[nmdbtag[i]+'T'].max()
            # if df[nmdbtag[i]+'T'].min()<ymin: ymin=df[nmdbtag[i]+'T'].min()
         for i in range(N,Nall):
            if(0 < dfCur[nmdbtag[i]+'T'].count()): 
               axesT.plot(dfCur.index.values,Fact[i]*dfCur[nmdbtag[i]+'T'],'-',linewidth=0.8,label='{0:s} {1:s}'.format(Labels[i],sFact[i]))
               


         # plt.tick_params(axis='x', which='major', labelsize=0,direction='in',length=6)
         # plt.tick_params(axis='y', which='major', labelsize=0,direction='in',length=6)
         # plt.tick_params(axis='y', which='minor', labelsize=0,direction='in',length=3)

         # axesT.xaxis.set_major_locator(mdates.HourLocator(interval=xTickMajorHours))
         # axesT.xaxis.set_minor_locator(mdates.MinuteLocator(byminute=range(0,xTickMajorHours*60,xTickMajorHours*15)))
         # axesT.xaxis.set_major_formatter(mdates.DateFormatter('%m/%d\n%H:%M'))
         # axesT.set_xlim(now - timedelta(hours=12),now)
         # axesT.set_xlim(startTime,endTime)
         # if (0==yminT):
         #    axesT.set_ylim(yminT,(1.+limMargin)*ymaxT)
         # else:
         #    axesT.set_ymargin(limMargin)
         #    axesT.set_ylim(yminT,ymaxT)
         axesT.set_ylim(yminT,ymaxT)
         # axesT.set_ylim(yminT,ymaxT-1)
         axesT.set_ylabel('Rate [count / minute]\n3-min moving average',fontsize=fontsize+1)
         plt.tick_params(axis='x', which='major', labelsize=0,direction='in',length=6)
         plt.tick_params(axis='x', which='minor', labelsize=0,direction='in',length=3)
         plt.tick_params(axis='y', which='major', labelsize=fontsize,direction='in',length=6)

         plt.grid(axis='both',which='both',linewidth=0.5,linestyle=':',color='gray')

         plt.legend(bbox_to_anchor=(1.01,0.525), loc="center left", borderaxespad=0,
                  fontsize=fontsize,labelspacing=0.0,frameon=False)
         # plt.title('GLE Alarm from Bartol Neutron Monitors - {0:s} UT'.format(dfCur.index.values[-1]),fontsize=fontsize)
         # plt.title('GLE Alarm from Bartol Neutron Monitors - {0:s} {1:s} {2:s} UT'.format(df.index.values[r],r'$@$',r.strftime("%H:%M:%S")),fontsize=fontsize)


      # axesal = fig.add_subplot(5,1,(1,1),sharex=axes)
      # axesal = fig.add_subplot(pAll,1,(pAll-4,pAll-4),sharex=axesT)
      axesal = fig.add_subplot(pAll,1,(pAll-(2+pT),pAll-(2+pT)),sharex=axes)
      dfStatus = dfCur[dfCur['Status']==3]
      # print(dfStatus) #DEBUG
      # plt.tick_params(axis='x', which='major', labelsize=0,direction='in',length=6)
      # plt.tick_params(axis='x', which='minor', labelsize=0,direction='in',length=3)
      # plt.tick_params(axis='y', which='major', labelsize=0,direction='in',length=0)
      # plt.tick_params(axis='y', which='major', labelsize=0,direction='in',length=6)
      if networkAwareAlert :
         # axesal.bar(df.index.values, df['Bartol_Above'], label='Bartol Simpson')
         # axesal.plot(df.index.values, df['Bartol_Above'], label='Bartol Simpson')
         axesal.set_ylim(0,ymaxAl+0.75)
         axesal.fill_between(x=df.index.values, y1=0, y2=0.5, color='lightgrey', alpha=0.2)
         axesal.fill_between(x=df.index.values, y1=0.5, y2=1.5, color='lightblue', alpha=0.2)
         axesal.fill_between(x=df.index.values, y1=1.5, y2=2.5, color='lightyellow', alpha=0.2)
         axesal.fill_between(x=df.index.values, y1=2.5, y2=ymaxAl+0.75, color='pink', alpha=0.2)
         # axesal.stackplot(dfCur.index.values, [dfCur['Bartol_Above'],dfCur['Extended_Above'],dfCur['Intl_Above']] , step="post", labels=['Current', 'Full Simpson', 'International'])
         axesal.stackplot(dfCur.index.values, [dfCur['Bartol_Above'],dfCur['Extended_Above'],dfCur['Intl_Above']] , step="post", labels=['Prototype', '+ Simpson Network', '+ Mawson'])
         # axesal.stackplot(dfCur.index.values, [dfCur['Bartol_Above'],dfCur['Extended_Above']] , labels=['Bartol Simpson', 'Extended Simpson'])
         # axesal.bar(dfCur.index.values, dfCur['Bartol_Above'],label='Bartol Simpson')
         # axesal.bar(df.index.values, df['Bartol_Above'], bottom=df['Bartol_Above'], label='Extended Simpson')
         # axesal.bar(df.index.values, df['Bartol_Above'], bottom=df['Bartol_Above']+df['Extended_Above'], label='International')
      else :
         
         axesal.plot(dfStatus.index.values,3*np.ones(len(dfStatus)),'o',color='red',label=Status[3])
         dfStatus = dfCur[dfCur['Status']==2]
         # plt.tick_params(axis='x', which='major', labelsize=0,direction='in',length=6)
         # plt.tick_params(axis='x', which='minor', labelsize=0,direction='in',length=3)
         # plt.tick_params(axis='y', which='major', labelsize=0,direction='in',length=6)
         axesal.plot(dfStatus.index.values,2*np.ones(len(dfStatus)),'o',color='orange',label=Status[2])
         dfStatus = dfCur[dfCur['Status']==1]
         # plt.tick_params(axis='x', which='major', labelsize=0,direction='in',length=6)
         # plt.tick_params(axis='x', which='minor', labelsize=0,direction='in',length=3)
         # plt.tick_params(axis='y', which='major', labelsize=0,direction='in',length=6)
         axesal.plot(dfStatus.index.values,1*np.ones(len(dfStatus)),'o',color='blue',label=Status[1])
         dfStatus = dfCur[dfCur['Status']==0]
         # plt.tick_params(axis='x', which='major', labelsize=0,direction='in',length=6)
         # plt.tick_params(axis='x', which='minor', labelsize=0,direction='in',length=3)
         # plt.tick_params(axis='y', which='major', labelsize=0,direction='in',length=6)
         axesal.plot(dfStatus.index.values,0*np.ones(len(dfStatus)),'o',color='gray',label=Status[0])
         axesal.set_ylim(0,3.75)
         # axesal.set_ylabel('Alarm Level',fontsize=fontsize+1)
         axesal.set_ylabel('Alarm Level',fontsize=fontsize)

      plt.tick_params(axis='x', which='major', labelsize=0,direction='in',length=6)
      plt.tick_params(axis='x', which='minor', labelsize=0,direction='in',length=3)
      # plt.tick_params(axis='y', which='major', labelsize=fontsize,direction='in',length=6)
      # plt.tick_params(axis='y', which='minor', labelsize=0,direction='in',length=3)
      # axesal.plot(df[df['Status']==3]['Time'],3*np.ones(len(df[df['Status']==3])),'o',color='red',label=Status[3])
      # axesal.plot(df[df['Status']==2]['Time'],2*np.ones(len(df[df['Status']==2])),'o',color='orange',label=Status[2])
      # axesal.plot(df[df['Status']==1]['Time'],1.*np.ones(len(df[df['Status']==1])),'o',color='blue',label=Status[1])
      # axesal.plot(df[df['Status']==0]['Time'],0*np.ones(len(df[df['Status']==0])),'o',color='gray',label=Status[0])
      # axesal.set_ylabel('Alarm Level',fontsize=fontsize+1)
      plt.legend(bbox_to_anchor=(1.01,0.48), loc="center left", borderaxespad=0,
               fontsize=fontsize,labelspacing=0.5,frameon=False)

      # print(axes.get_xticklabels())  #DEBUG
      if lenGP > 0:
         # axesGP = fig.add_subplot(pAll,1,(pAll-5,pAll-5),sharex=axesT)
         axesGP = fig.add_subplot(pAll,1,(pAll-(3+pT),pAll-(3+pT)),sharex=axes)
         axesGP.set_ylim(yminGP,ymaxGP)
         # axesGP.yaxis.set_minor_locator(subs='auto')
         axesGP.set_yscale('log')
         # axesGP.yaxis.set_minor_locator(mticker.LogLocator( numticks=10, subs='auto'))
         # print(axesGP.yaxis.get_tick_space()) #DEBUG
         llGP=mticker.LogLocator(base=10.0, numticks=int(math.ceil(ydeltaGP)))#subs=np.arange(2, 10) * 0.1)
         llmGP=mticker.LogLocator(base=10.0, numticks=int(math.ceil(ydeltaGP))*9, subs='auto')#subs=np.arange(2, 10) * 0.1)
         axesGP.yaxis.set_major_locator(llGP)
         axesGP.yaxis.set_minor_locator(llmGP)
         # axesGP.yaxis.set_minor_formatter(mticker.LogFormatterMathtext(base=10.0,  labelOnlyBase=False))
         axesGP.yaxis.set_major_formatter(mticker.LogFormatterMathtext(base=10.0,  labelOnlyBase=False,minor_thresholds=(0, 0)))
         # print(axesGP.get_yticks()) #DEBUG
         axesGP.plot(dfGPCur.index.values,dfGPCur['p3_flux_ic'],color='darkred',label='>=10 MeV')
         axesGP.plot(dfGPCur.index.values,dfGPCur['p7_flux_ic'],color='darkblue',label='>=100 MeV')
         # axesGP.plot(df500['time_tag2'].to_numpy(),df500['flux'].to_numpy(),color='pink',label='>=500 MeV')
         # axesGP.set_ymargin(10*limMargin)
         # axesGP.yaxis.set_major_locator(mticker.LogLocator(  subs='auto'))
         axesGP.set_ylabel('GOES\n Particles\n cm$^{-2}$s$^{-1}$sr$^{-1}$ ',fontsize=fontsize,color='k', multialignment='center')
         plt.tick_params(axis='x', which='major', labelsize=0,direction='in',length=6)
         plt.tick_params(axis='x', which='minor', labelsize=0,direction='in',length=3)
         plt.tick_params(axis='y', which='major', labelsize=fontsize,direction='in',length=6)
         plt.tick_params(axis='y', which='minor', labelsize=0,direction='in',length=3)
         # print(axesGP.yaxis.get_tick_params(which='minor')) #DEBUG
         plt.grid(axis='x',which='major',linewidth=0.5,linestyle='-',color='gray')
         plt.grid(axis='x',which='minor',linewidth=0.5,linestyle=':',color='gray')
         plt.grid(axis='y',which='major',linewidth=0.5,linestyle=':',color='gray')
         # plt.legend(loc='upper right',fontsize=fontsize,ncol=1)
         plt.legend(bbox_to_anchor=(1.01,0.48), loc="center left", borderaxespad=0,
               fontsize=fontsize,labelspacing=0.5,frameon=False)
         # print(ymaxGP) #DEBUG
         # axesGP.set_xticklabels([])
         if alarmLineGPShow :
   
            alarmLineGP=datetime.combine(startDay, datetime.min.time())+timedelta(hours=10,minutes=29)
            if (dfCur.index[-1]>=alarmLineGP):
               # print(alarmLineGP) #DEBUG

               axesGP.axvline(datetime.combine(startDay, datetime.min.time())+timedelta(hours=10,minutes=29),color=alarmColors[2])
      if lenGX > 0:
            # axes = fig.add_subplot(23,1,(1,2))
            # axes.plot(dfx['time_tag2'].to_numpy(),dfx['flux'].to_numpy(),'-',color='k')
            # plt.tick_params(axis='x', which='major', labelsize=0,direction='in',length=6)
            # plt.tick_params(axis='x', which='minor', labelsize=0,direction='in',length=3)
            # plt.tick_params(axis='y', which='major', labelsize=fontsize,direction='in',length=6)
            # plt.tick_params(axis='y', which='minor', labelsize=0,direction='in',length=3)
            # plt.grid(axis='x',which='major',linewidth=0.5,linestyle='-',color='gray')
            # plt.grid(axis='x',which='minor',linewidth=0.5,linestyle=':',color='gray')
            # plt.grid(axis='y',which='major',linewidth=0.5,linestyle=':',color='gray')
            # axes.set_yscale('log')
            # axes.set_ylabel('GOES\n X-ray flux\n Watts m$^{-2}$ ',fontsize=fontsize,color='k', multialignment='center')
            # axes.text(enddisp,1e-4,' X',fontsize=fontsize,color='k',verticalalignment='center', horizontalalignment='left')
            # axes.text(enddisp,1e-5,' M',fontsize=fontsize,color='k',verticalalignment='center', horizontalalignment='left')
            # axes.text(enddisp,1e-6,' C',fontsize=fontsize,color='k',verticalalignment='center', horizontalalignment='left')
            # axes.set_xlim(startdisp,enddisp)
   
         # axesGX = fig.add_subplot(pAll,1,(pAll-(4+pG),pAll-(4+pG)),sharex=axesT)
         axesGX = fig.add_subplot(pAll,1,(pAll-(2+pT+pG),pAll-(2+pT+pG)),sharex=axes)
         axesGX.plot(dfGXCur.index.values,dfGXCur['xs'],color='green',label='XS')
         axesGX.plot(dfGXCur.index.values,dfGXCur['xl'],color='k',label='XL')
         # axesGP.plot(df500['time_tag2'].to_numpy(),df500['flux'].to_numpy(),color='pink',label='>=500 MeV')
         # axesGX.set_ymargin(10*limMargin)
         axesGX.set_ylim(yminGX,ymaxGX)
         axesGX.set_yscale('log')
         axesGX.set_ylabel('GOES\n X-ray flux\n Watts m$^{-2}$ ',fontsize=fontsize,color='k', multialignment='center')
         plt.tick_params(axis='x', which='major', labelsize=0,direction='in',length=6)
         plt.tick_params(axis='x', which='minor', labelsize=0,direction='in',length=3)
         plt.tick_params(axis='y', which='major', labelsize=fontsize,direction='in',length=6)
         plt.tick_params(axis='y', which='minor', labelsize=0,direction='in',length=3)
         plt.grid(axis='x',which='major',linewidth=0.5,linestyle='-',color='gray')
         plt.grid(axis='x',which='minor',linewidth=0.5,linestyle=':',color='gray')
         plt.grid(axis='y',which='major',linewidth=0.5,linestyle=':',color='gray')
         # plt.legend(loc='upper right',fontsize=fontsize,ncol=1)
         plt.legend(bbox_to_anchor=(1.01,0.48), loc="center left", borderaxespad=0,
               fontsize=fontsize,labelspacing=0.5,frameon=False)
         # print(pAll) #DEBUG
         # print(axesGX.get_yticks()) #DEBUG
         # print(axesGX.get_yticklabels()) #DEBUG
         # sys.exit() #DEBUG

         # axesGX.set_xticklabels([])

      # axes.axhline(y=4,color='red')
      axes.fill_between(x=df.index.values, y1=yminI, y2=4.0, color='lightgrey', alpha=0.5)

      for i in range(len(alarmLines)):
         if dfCur.index.values[-1] >= alarmLines[i]:
            axes.axvline(alarmLines[i],color=alarmColors[i])
            if (ratePlot): axesT.axvline(alarmLines[i],color=alarmColors[i])
      if showBaselines :
         for i in range(len(baselines)):
            axes.axvline(baselines[i],color='green')
            if (ratePlot):axesT.axvline(baselines[i],color='green')
         if (3 > dfCur.iloc[-1]['Status']) :
            baselines = [df.index[r-85],df.index[r-10] ]

      # print(axes.get_xticklabels())  #DEBUG
      # if (ratePlot):axesT.set_xticklabels([])

      # axesal.set_xticklabels([])
      # axesal.set_yticks([])
      # axesal.set_yticklabels([])
      # plt.yticks(np.arange(0, 4, 1.0))

      # axesal.xaxis.set_major_locator(mdates.HourLocator(interval=1))
      # axesal.xaxis.set_major_formatter(mdates.DateFormatter('%m/%d\n%H:%M'))
      plt.grid(axis='x',which='both',linewidth=0.5,linestyle=':',color='gray')
      # plt.legend(bbox_to_anchor=(1.01,0.48), loc="center left", borderaxespad=0,
      #          fontsize=fontsize,labelspacing=0.5,frameon=False)

      # axesal.text(axesal.get_xlim()[0] + 0.02*(axesal.get_xlim()[1] -axesal.get_xlim()[0] ) ,
      #          axesal.get_ylim()[1]- 0.15*(axesal.get_ylim()[1] -axesal.get_ylim()[0] ),
      #          # axesal.get_ylim()[1]- 0.12*(axesal.get_ylim()[1] -axesal.get_ylim()[0] ),
      #          "Current: ", horizontalalignment='left', fontsize=fontsize+1,zorder=10)

      # axesal.text(axesal.get_xlim()[0] + 0.11*(axesal.get_xlim()[1] -axesal.get_xlim()[0] ) ,
      #          axesal.get_ylim()[1]- 0.15*(axesal.get_ylim()[1] -axesal.get_ylim()[0] ),
      #          # axesal.get_ylim()[1]- 0.12*(axesal.get_ylim()[1] -axesal.get_ylim()[0] ),
      #          "{0:s}".format(Status[int(LastStatus)]), horizontalalignment='left',color=Statuscol[int(LastStatus)], fontsize=fontsize+1,zorder=10)

      axesal.set_ylabel("Number of stations\nabove threshold",fontsize=fontsize,multialignment='center')

      # axesal.text(axesal.get_xlim()[0], axesal.get_ylim()[1]+ 0.017*(axesal.get_ylim()[1] -axesal.get_ylim()[0] ),
            # 'Last update: {0:s} UT'.format(dfCur.index.values[-1]),fontsize=fontsize+1,horizontalalignment='left')


      plt.subplots_adjust(left=0.1, bottom=0.06, right=0.8, top=0.95, wspace=0, hspace=0.00)

      # print(yminT,ymaxT,yminI,ymaxI,yminGP,ymaxGP,yminGX,ymaxGX) #DEBUG
      # print(axesT.yaxis.get_data_interval(),axes.yaxis.get_data_interval(),axesGP.yaxis.get_data_interval(),axesGX.yaxis.get_data_interval()) #DEBUG
      # axesLabels = axes.get_xticklabels()
      # for lbl in axesLabels :
      #    lbl.set_visible(True)
      # print(axes.get_xticklabels())


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
