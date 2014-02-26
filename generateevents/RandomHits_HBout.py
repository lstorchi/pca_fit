#!/usr/bin/python

import sys, string, os, time, re
import random

def main():
  print "Generating random hits for gf test"
  if len(sys.argv)<5:
    print "usage: RandomHits.py numEvts numRoads numHitsPerLayer outFile"
    return -1
 
  #the following two numbers are fixed		
  NSVXhits = 5
  NSVThits = 6
 
  Ntest = int(sys.argv[1])	
  NRoads = int(sys.argv[2])	
  NHitsPerLayer = int(sys.argv[3])	
  fout = open(sys.argv[4],"w")
  fout2 = open("hitsDecoded","w")
  for ntest in range(0,Ntest):
    hitList = list()
    hitListDecoded = list()
    print "***** Event  #  ",
    print ntest
    for h in range(0,NRoads): 
      print "Road   #  ",
      print NRoads
      for i in range(0,NSVThits): 
        print "hit num #  ",  
        print i
        #for each layer I generate a random number of hits
        #for j in range(random.randint(1,NHitsPerLayer)):
        #for each layer I generate a fixed number of hits
        for j in range(0,NHitsPerLayer):
          print "  same layer hit num #  ",  
          print j
          if i < NSVXhits:
            # 15 bits for hit coordinate (bit 0 - 14)
            randNum = random.randint(0, 0x7fff)
            hit = randNum
            print >> fout2,"hit %.6x" % hit  
            #hit = 1
            #print "%.6x" % hit
            #3 bits for the barrel (random number) bit 15-17 (see pg 73 rafanelli) 
            zbarrel = random.randint(0,0x3)
            print >> fout2,"zbarrel %.6x" % zbarrel  
            #print "zbarrel: %x" % zbarrel
            hit +=  (zbarrel << 15)  
            #print "%.6x" % hit
            # 3 bits for the layer, bit 18-20
            layer = i
            print >> fout2,"layer %.6x" % layer  
            hit +=  (layer << 18)
            #print "%.6x" % hit
            if layer != 0:
              print >> fout,"%.6x" % hit
            hitList.append(hit)	
          
          elif i == NSVXhits:
            # XFT information for the AM, 18 bits (0-17)
            hit = 0x1feed
            print >> fout2,"XFT for AM %.6x" % hit  
            #print "%.6x" % hit	
            # layer (XFT), 3 bits, (18-20)
            layer = i
            print >> fout2,"layer(XFT) %.6x" % layer
            hit += (layer << 18)
            #print "%.6x" % hit
            print >> fout,"%.6x" % hit
            hitList.append(hit)
            
            #phi, 12 bits (0-11)
            phi = random.randint(0,0xfff)
            hit = phi
            print >> fout2,"phi %.6x" % phi
            #print "%.6x" % hit
            #curvature, 7 bits (last bit is the sign), 12-18
            curv = random.randint(0,0x3f)
            print >> fout2,"curv %.6x" % curv
            hit += (curv << 12)
            #print "%.6x" % hit
            print >> fout,"%.6x" % hit
            hitList.append(hit)
     
      #AM road ID, 21 bits, 0-20
      amroadID = random.randint(0,0x1fffff)
      hit = amroadID
      print >> fout2,"amroadID %.6x" % amroadID
      # print "%.6x" % hit
      #EP = 1
      hit += (1<<21)
      print >> fout2,"EP %.6x" % 1
      #print "%.6x" % hit
      print >> fout,"%.6x" % hit
      hitList.append(hit)
          
    #last word (must contain EE = 1)
    #for the moment,21 bits random
    lastword = random.randint(0,0x1fffff)
    hit = lastword
    print >> fout2,"lastword %.6x" % lastword
    #print "%.6x" % hit
    #EP = 1 (bit 21)
    hit += (1<<21)
    print >> fout2,"EP %.6x" % 1
    #print "%.6x" % hit
    #EE = 1 (bit 22)
    hit += (1<<22)
    #print "%.6x" % hit
    print >> fout2,"EE %.6x" % 1
    print >> fout,"%.6x" % hit
    hitList.append(hit)
  
  #fout.write("\n")
  hitList.append(hit)

if __name__=="__main__":
  main()
