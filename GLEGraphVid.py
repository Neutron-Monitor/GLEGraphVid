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
# 1.1.0 Add GOES data
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
   urlalarm='./GLE_Alarm.png' #DEBUG
   initMinutes=30
   frameNum = 0
   fileGOESProton =''
   fileGOESXray =''

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

   try:
      opts, args = getopt.getopt(argv,"hr:s:e:i:o:p:x:")
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
   Labels=  ['Inuvik','Fort Smith','Peawanuck','Nain','Newark','Thule','South Pole','South Pole - bare','McMurdo','Jang Bogo']
   InAlert= [1       ,1           ,1          ,1     ,0       ,1      ,1                       ,0                  ,1        ,0          ]
   sFact= ['','',' *2','',' *2','',' /2','','','']
   Fact=  [1.,1.,2.   ,1.,2.    ,1.,0.5 ,1.,1.,1.]
#TODO remove not inalert
#TODO calc vert lines for events before hand show later

   #History from Makejson_ql.py
   #History= [0.597135,(0.598/0.94696),1.35333,0.59686,0.54518,0.57732*0.6,0.705,0.52308,0.52308]
   History= [0.597135,0.598,1.35333,0.59686,0.54518,0.57732*0.6,0.52308,0.52308,0.705]
   History=np.array(History)
   History=History/0.6

   lenGP=0
   lenGX=0

   df = pd.read_csv('{0:s}/GLE_Day_{1:s}.txt'.format(
                     Outpath,startDay.strftime("%Y%m%d")),
                     sep=',',date_format='%y/%m/%d %H:%M:%S',
                     index_col=0)
   df=df.drop(columns=['Time.1'])
   # print(df)  #DEBUG
   # print(df.info(verbose=True, show_counts=True))  #DEBUG


   #Number of stations
   N= len(nm)
   Notused='$^{\dagger}$Not used'
   i=1
   while i < N:
      if not InAlert[i]:
         nm.append(nm.pop(i)    )
         nmdbtag.append(nmdbtag.pop(i))
         Labels.append('$^{\dagger}'+Labels.pop(i) )
         InAlert.append(InAlert.pop(i))
         sFact.append(sFact.pop(i))
         Fact.append(Fact.pop(i))
         N-=1
         #TODO add df analysis for no alert
      else : i+=1
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
   # print(df)  #DEBUG
   # print(df.info(verbose=True, show_counts=True))  #DEBUG
   # sys.exit() #DEBUG

   lineGP=454
   if os.path.isfile(fileGOESProton):
      dfGP = pd.read_csv(fileGOESProton,
                     sep=',',date_format='%y-%m-%d %H:%M:%S.%f',parse_dates=['time_tag'],
                     index_col=0, skiprows=lineGP, na_values=np.nan)
                     # usecols=[26,30])
                     # usecols=['p3_flux_ic','p7_flux_ic'])
      dfGP.index = pd.to_datetime(dfGP.index)

      for c in ['e1_flux_ic','e2_flux_ic','e3_flux_ic','p1_flux','p2_flux','p3_flux','p4_flux','p5_flux','p6_flux','p7_flux','a1_flux','a2_flux','a3_flux','a4_flux','a5_flux','a6_flux','p1_flux_c','p2_flux_c','p3_flux_c','p4_flux_c','p5_flux_c','p6_flux_c','p7_flux_c','p1_flux_ic','p2_flux_ic','p4_flux_ic','p5_flux_ic','p6_flux_ic']:
         dfGP=dfGP.drop(columns=c) #drop all but p3_flux_ic, p7_flux_ic
      lenGP=len(dfGP)
      # print(dfGP)  #DEBUG
      # print(dfGP.info(verbose=True, show_counts=True))  #DEBUG
   # sys.exit() #DEBUG

   lineX=118
   if os.path.isfile(fileGOESXray):
      dfGX = pd.read_csv(fileGOESXray,
                     sep=',',date_format='%y-%m-%d %H:%M:%S.%f',parse_dates=['time_tag'],
                     index_col=0, skiprows=lineX, na_values=np.nan)
                     # usecols=[26,30])
                     # usecols=['p3_flux_ic','p7_flux_ic'])
      dfGX.index = pd.to_datetime(dfGX.index)

      # for c in ['e1_flux_ic','e2_flux_ic','e3_flux_ic','p1_flux','p2_flux','p3_flux','p4_flux','p5_flux','p6_flux','p7_flux','a1_flux','a2_flux','a3_flux','a4_flux','a5_flux','a6_flux','p1_flux_c','p2_flux_c','p3_flux_c','p4_flux_c','p5_flux_c','p6_flux_c','p7_flux_c','p1_flux_ic','p2_flux_ic','p4_flux_ic','p5_flux_ic','p6_flux_ic']:
      #    dfGP=dfGP.drop(columns=c) #drop all but p3_flux_ic, p7_flux_ic
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
   pG=0
   if lenGP>0:
      ymaxGP=200.
      yminGP=0.1
      pG=1
      # print(dfGP['p3_flux_ic'].max())
      if dfGP['p3_flux_ic'].max()>ymaxGP: ymaxGP=dfGP['p3_flux_ic'].max()
      # if dfGP['p3_flux_ic'].min()<yminGP: yminGP=dfGP'p3_flux_ic'].min()
      if dfGP['p7_flux_ic'].max()>ymaxGP: ymaxGP=dfGP['p7_flux_ic'].max()
      # if dfGP['p7_flux_ic'].min()<yminGP: yminGP=dfGP['p7_flux_ic'].min()


   if lenGX>0:
      ymaxGX=1e-9
      yminGX=1e-11
      if dfGX['xs'].max()>ymaxGX: ymaxGX=dfGX['xs'].max()
      # if dfGX['xs'].min()<yminGX: yminGX=dfGX'xs'].min()
      if dfGX['xl'].max()>ymaxGX: ymaxGX=dfGX['xl'].max()
      # if dfGP['xl'].min()<yminGP: yminGP=dfGP['xl'].min()

      pG+=1
   pAll=5+pG
   LastStatus=0
   for r in range(initMinutes, endMinutes-startMinutes) :
   # for r in range(endMinutes-startMinutes-1, endMinutes-startMinutes) :

      dfCur=df.iloc[0:r]

      if lenGP>0:
         dfGPCur = dfGP[startTime:dfCur.index.values[-1]]
         # print(dfGPCur)  #DEBUG
         # print(dfGPCur.info(verbose=True, show_counts=True))  #DEBUG
      if lenGX>0:
         dfGXCur = dfGX[startTime:dfCur.index.values[-1]]
         # print(dfGXCur)  #DEBUG
         # print(dfGXCur.info(verbose=True, show_counts=True))  #DEBUG

      fig=plt.figure(figsize=(14, 11), dpi=80)

      # axesT = fig.add_subplot(5,1,(2,3))
      axesT = fig.add_subplot(pAll,1,(pAll-3,pAll-2))

      for i in range(N):
         # axesT.plot(df['Time'],Fact[i]*df[nmdbtag[i]+'T'],'-',linewidth=0.8,label=r'{0:s} {1:s}'.format(Labels[i],sFact[i]))
         axesT.plot(dfCur.index.values,Fact[i]*dfCur[nmdbtag[i]+'T'],'-',linewidth=0.8,label='{0:s} {1:s}'.format(Labels[i],sFact[i]))
         # if df[nmdbtag[i]+'T'].max()>ymax: ymax=df[nmdbtag[i]+'T'].max()
         # if df[nmdbtag[i]+'T'].min()<ymin: ymin=df[nmdbtag[i]+'T'].min()

      # plt.tick_params(axis='x', which='major', labelsize=0,direction='in',length=6)
      # plt.tick_params(axis='y', which='major', labelsize=fontsize,direction='in',length=6)
      # plt.tick_params(axis='y', which='major', labelsize=fontsize,direction='in',length=6)

      # axesT.xaxis.set_major_locator(mdates.HourLocator(interval=1))
      # axesT.xaxis.set_minor_locator(mdates.MinuteLocator(interval=15))
      # axesT.xaxis.set_major_formatter(mdates.DateFormatter('%m/%d\n%H:%M'))
      # axesT.set_xlim(now - timedelta(hours=12),now)
      axesT.set_xlim(startTime,endTime)
      axesT.set_ylim(yminT,ymaxT-1)
      axesT.set_ylabel('Rate [count / minute]\n3-min moving average',fontsize=fontsize+1)

      # plt.grid(axis='both',which='both',linewidth=0.5,linestyle=':',color='gray')

      plt.legend(bbox_to_anchor=(1.01,0.525), loc="center left", borderaxespad=0,
               fontsize=fontsize,labelspacing=1,frameon=False)
      # plt.title('GLE Alarm from Bartol Neutron Monitors - {0:s} UT'.format(dfCur.index.values[-1]),fontsize=fontsize)
      # plt.title('GLE Alarm from Bartol Neutron Monitors - {0:s} {1:s} {2:s} UT'.format(df.index.values[r],r'$@$',r.strftime("%H:%M:%S")),fontsize=fontsize)

      # axes = fig.add_subplot(5,1,(4,5),sharex=axesT)
      axes = fig.add_subplot(pAll,1,(pAll-1,pAll))


      for i in range(N):
         axes.plot(dfCur.index.values,100.*(dfCur[nmdbtag[i]+'Ith']-1.),'-',linewidth=0.8,label='{0:s}'.format(Labels[i]))


      # plt.tick_params(axis='x', which='major', labelsize=fontsize+1,direction='in',length=6)
      # plt.tick_params(axis='y', which='major', labelsize=fontsize,direction='in',length=6)
      # plt.tick_params(axis='y', which='major', labelsize=fontsize,direction='in',length=6)

      axes.xaxis.set_major_locator(mdates.HourLocator(interval=1))
      axes.xaxis.set_minor_locator(mdates.MinuteLocator(interval=15))
      axes.xaxis.set_major_formatter(mdates.DateFormatter('%m/%d\n%H:%M'))
      ###axes.set_xlim(now - timedelta(hours=8),now)
      axes.set_xlim(startTime,endTime)
      axes.set_ylim(yminI,ymaxI-0.001)
      axes.set_ylabel('Rate increase [%]\n(3-min tr. moving average)',fontsize=fontsize+1)
      axes.axhline(y=dfCur.iloc[-1]['Status'],linewidth=0.5,linestyle='-',color='red')
      # axes.axhline(y=Level,linewidth=0.5,linestyle='-',color='red')

      # plt.grid(axis='both',which='major',linewidth=0.5,linestyle=':',color='gray')
      plt.legend(bbox_to_anchor=(1.01,0.525), loc="center left", borderaxespad=0,
               fontsize=fontsize,labelspacing=1,frameon=False)
      #plt.title('GLE Alarm from Bartol Neutron Monitors - Last update: {0:s} {1:s} {2:s} UT'.format(now.strftime("%Y-%m-%d"),r'$@$',now.strftime("%H:%M:%S")),fontsize=fontsize)
      # plt.title('GLE Alarm from Bartol Neutron Monitors - Last update: {0:s} UT'.format(dfCur.index.values[-1]),fontsize=fontsize)

      axes.text(axes.get_xlim()[1] + 0.19*(axes.get_xlim()[1] -axes.get_xlim()[0] ) ,
               axes.get_ylim()[0]+ 0.0*(axes.get_ylim()[1] -axes.get_ylim()[0] ),
               Notused, horizontalalignment='left', fontsize=fontsize-2,zorder=10)


      # axesal = fig.add_subplot(5,1,(1,1),sharex=axes)
      axesal = fig.add_subplot(pAll,1,(pAll-4,pAll-4),sharex=axesT)
      dfStatus = dfCur[df['Status']==3]
      # print(dfStatus) #DEBUG
      axesal.plot(dfStatus.index.values,3*np.ones(len(dfStatus)),'o',color='red',label=Status[3])
      dfStatus = dfCur[df['Status']==2]
      axesal.plot(dfStatus.index.values,2*np.ones(len(dfStatus)),'o',color='orange',label=Status[2])
      dfStatus = dfCur[df['Status']==1]
      axesal.plot(dfStatus.index.values,1*np.ones(len(dfStatus)),'o',color='blue',label=Status[1])
      dfStatus = dfCur[df['Status']==0]
      axesal.plot(dfStatus.index.values,0*np.ones(len(dfStatus)),'o',color='gray',label=Status[0])
      axesal.set_ylim(0,3.75)
      axesal.set_ylabel('Alarm Level',fontsize=fontsize+1)
      # axesal.plot(df[df['Status']==3]['Time'],3*np.ones(len(df[df['Status']==3])),'o',color='red',label=Status[3])
      # axesal.plot(df[df['Status']==2]['Time'],2*np.ones(len(df[df['Status']==2])),'o',color='orange',label=Status[2])
      # axesal.plot(df[df['Status']==1]['Time'],1.*np.ones(len(df[df['Status']==1])),'o',color='blue',label=Status[1])
      # axesal.plot(df[df['Status']==0]['Time'],0*np.ones(len(df[df['Status']==0])),'o',color='gray',label=Status[0])
      # axesal.set_ylabel('Alarm Level',fontsize=fontsize+1)

      if lenGP > 0:
         axesGP = fig.add_subplot(pAll,1,(pAll-5,pAll-5),sharex=axesT)
         axesGP.plot(dfGPCur.index.values,dfGPCur['p3_flux_ic'],color='darkred',label='>=10 MeV')
         axesGP.plot(dfGPCur.index.values,dfGPCur['p7_flux_ic'],color='darkblue',label='>=100 MeV')
         # axesGP.plot(df500['time_tag2'].to_numpy(),df500['flux'].to_numpy(),color='pink',label='>=500 MeV')
         # plt.tick_params(axis='x', which='major', labelsize=0,direction='in',length=6)
         # plt.tick_params(axis='x', which='minor', labelsize=0,direction='in',length=3)
         # plt.tick_params(axis='y', which='major', labelsize=fontsize,direction='in',length=6)
         # plt.tick_params(axis='y', which='minor', labelsize=0,direction='in',length=3)
         # plt.grid(axis='x',which='major',linewidth=0.5,linestyle='-',color='gray')
         # plt.grid(axis='x',which='minor',linewidth=0.5,linestyle=':',color='gray')
         # plt.grid(axis='y',which='major',linewidth=0.5,linestyle=':',color='gray')
         axesGP.set_ylim(yminGP,ymaxGP)
         axesGP.set_yscale('log')
         axesGP.set_ylabel('GOES\n Particles\n cm$^{-2}$s$^{-1}$sr$^{-1}$ ',fontsize=fontsize,color='k', multialignment='center')
         # plt.legend(loc='upper right',fontsize=fontsize,ncol=1)
         # print(ymaxGP) #DEBUG
         axesGP.set_xticklabels([])
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
   
         axesGX = fig.add_subplot(pAll,1,(pAll-(4+pG),pAll-(4+pG)),sharex=axesT)
         axesGX.plot(dfGXCur.index.values,dfGXCur['xs'],color='k',label='XS')
         axesGX.plot(dfGXCur.index.values,dfGXCur['xl'],color='orange',label='XL')
         # axesGP.plot(df500['time_tag2'].to_numpy(),df500['flux'].to_numpy(),color='pink',label='>=500 MeV')
         # plt.tick_params(axis='x', which='major', labelsize=0,direction='in',length=6)
         # plt.tick_params(axis='x', which='minor', labelsize=0,direction='in',length=3)
         # plt.tick_params(axis='y', which='major', labelsize=fontsize,direction='in',length=6)
         # plt.tick_params(axis='y', which='minor', labelsize=0,direction='in',length=3)
         # plt.grid(axis='x',which='major',linewidth=0.5,linestyle='-',color='gray')
         # plt.grid(axis='x',which='minor',linewidth=0.5,linestyle=':',color='gray')
         # plt.grid(axis='y',which='major',linewidth=0.5,linestyle=':',color='gray')
         axesGX.set_ylim(yminGX,ymaxGX)
         axesGX.set_yscale('log')
         axesGX.set_ylabel('GOES\n X-ray flux\n Watts m$^{-2}$ ',fontsize=fontsize,color='k', multialignment='center')
         # plt.legend(loc='upper right',fontsize=fontsize,ncol=1)
         # print(pAll) #DEBUG
         # sys.exit() #DEBUG

         axesGX.set_xticklabels([])

      plt.tick_params(axis='x', which='major', labelsize=0,direction='in',length=6)
      plt.tick_params(axis='x', which='minor', labelsize=0,direction='in',length=3)
      plt.tick_params(axis='y', which='major', labelsize=fontsize,direction='in',length=6)

      axesT.set_xticklabels([])
      axesal.set_xticklabels([])
      # plt.yticks(np.arange(0, 4, 1.0))

      # axesal.xaxis.set_major_locator(mdates.HourLocator(interval=1))
      # axesal.xaxis.set_major_formatter(mdates.DateFormatter('%m/%d\n%H:%M'))
      plt.grid(axis='x',which='both',linewidth=0.5,linestyle=':',color='gray')
      plt.legend(bbox_to_anchor=(1.01,0.48), loc="center left", borderaxespad=0,
               fontsize=fontsize,labelspacing=1,frameon=False)

      axesal.text(axesal.get_xlim()[0] + 0.02*(axesal.get_xlim()[1] -axesal.get_xlim()[0] ) ,
               axesal.get_ylim()[1]- 0.12*(axesal.get_ylim()[1] -axesal.get_ylim()[0] ),
               "Current: ", horizontalalignment='left', fontsize=fontsize+1,zorder=10)

      axesal.text(axesal.get_xlim()[0] + 0.11*(axesal.get_xlim()[1] -axesal.get_xlim()[0] ) ,
               axesal.get_ylim()[1]- 0.12*(axesal.get_ylim()[1] -axesal.get_ylim()[0] ),
               "{0:s}".format(Status[int(LastStatus)]), horizontalalignment='left',color=Statuscol[int(LastStatus)], fontsize=fontsize+1,zorder=10)


      # axesal.text(axesal.get_xlim()[0], axesal.get_ylim()[1]+ 0.017*(axesal.get_ylim()[1] -axesal.get_ylim()[0] ),
            # 'Last update: {0:s} UT'.format(dfCur.index.values[-1]),fontsize=fontsize+1,horizontalalignment='left')


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
