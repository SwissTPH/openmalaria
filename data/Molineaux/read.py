#!/usr/bin/env python3

# Purpose: read the list of asexual densities and calculate some derived numbers
# Author: Diggory Hardy
# Date: October 2013

# Input file is tab-separated with columns:
# Patient ID
# Day (starting 1, then every second day)
# Asexual parasite density in units of parasites per micro-litre

# Output file is semicolon-separated with columns:
# Parasite density at first local maximum, again parasites/μL
# Duration (last day with positive density minus first day with positive density)

class Reader:
    def __init__(self,out_stream):
        self.out_stream = out_stream
        self.reset_vars()
        self.patient=None
    
    def finish(self):
        self.flush_patient()
    
    def reset_vars(self):
        self.day = -1 # first day should be this plus 2 (i.e. 1)
        self.density = None
        self.first_patent = None
        self.last_patent = None
        self.first_local_max = None
    
    def flush_patient(self):
        if self.patient != None:
            assert self.first_patent == 1  # presumably
            assert self.last_patent > self.first_patent
            duration = self.last_patent - self.first_patent
            
            assert self.first_local_max > 0
            
            #print(self.patient,duration,self.first_local_max,sep="\t")
            print(duration,self.first_local_max,"// "+self.patient, sep=", ")
            print(duration, self.first_local_max, sep=";", file=self.out_stream)
    
    def read_line(self, line):
        row = line.split("\t")
        
        patient = row[0]
        if patient != self.patient:
            self.flush_patient()
            self.reset_vars()
            self.patient = patient
        
        day = int(row[1])
        assert day == self.day + 2
        self.day = day
        
        if row[2].strip() == '.':
            # no data; can't do much
            return
        
        density = int(row[2])
        
        if self.first_patent == None and density > 0:
            self.first_patent = day
        
        if density > 0:
            self.last_patent = day
        
        if self.first_local_max == None and self.density != None and self.density > density:
            self.first_local_max = self.density
        
        self.density = density

with open("35casesAsexDensities.txt") as in_stream, open("molineaux_params.csv", mode='w') as out:
    print("first local max parasites/microlitre","last day patent minus first", sep=";", file=out)
    
    in_stream.readline()  # discard header
    
    reader = Reader(out)
    for line in in_stream:
        reader.read_line(line)
    reader.finish()
