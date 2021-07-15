#!/usr/bin/python3.6
 
import sys
import os


def mapper():   
    filepath = os.environ["map_input_file"] 
    gfilename = os.path.split(filepath)[-1]  #get the names 
    for line in sys.stdin:
        line = line.strip()
        if line=="":
            continue
        temperature = line.split(',')[-1]
            # convert count (currently a string) to int
        print("%s\t%s" % (filename, temperature))

if __name__ == "__main__":
    mapper()