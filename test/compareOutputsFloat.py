#!/usr/bin/python
# -*- coding: utf-8 -*-
import sys
import string
import math

def main(*args):
  if (4!=len(args)):
    print "Usage: compareLogs.py logfile1 logfile2 n_prev_lines"
    return 1
  
  print "Comparing "+args[1]+" with "+args[2]
  file1=open(args[1], 'r')
  file2=open(args[2], 'r')
  pre_lines=int(args[3])+1
  line_count=0
  prev_lines1=[""]*pre_lines
  prev_lines2=[""]*pre_lines
  for line1 in file1:
    line_count+=1
    
    line_items1=string.split(line1)
    h_id1=int(line_items1[1])
    loc_id1=line_items1[2]
    value1=float(line_items1[3])
    
    line2=file2.readline()
    line_items2=string.split(line2)
    h_id2=int(line_items2[1])
    loc_id2=line_items2[2]
    value2=float(line_items2[3])
    
    if (loc_id1!=loc_id2) or (h_id1!=h_id2):
      print "Different summary outputs {0}:".format(line_count)
      print '-',line1,
      print '+',line2,
      return 2
    
    # Compare with relative precision. Although openMalaria normally rounds it's output to less decimal places than this...
    if math.fabs(value1 - value2) > 0.0000001 * max(math.fabs(value1),math.fabs(value2)):
      print "Different values output {0}".format(line_count)
      print '-',line1,
      print '+',line2,
      for j in range(0, pre_lines-1):
	print '-',prev_lines1[j],
	print '+',prev_lines2[j],
      return 1
    
    if(line_count % 100000 == 0):
      print (line_count)
    for i in range(pre_lines-1, 0, -1):
      prev_lines1[i]=prev_lines1[i-1]
      prev_lines2[i]=prev_lines2[i-1]
    prev_lines1[0]=line1
    prev_lines2[0]=line2
  print "No differences, ok..."
  return 0

if __name__ == '__main__':
  sys.exit(main(*sys.argv))
