#!/usr/bin/python3.6
 
from operator import itemgetter
import sys

current_name = None
current_temp = 0
filename = None
# input comes from STDIN
for line in sys.stdin:  
    # remove leading and trailing whitespace
    line = line.strip()  
    # parse the input we got from mapper.py
    filename, temp = line.split('\t', 1)  
    # convert temp (currently a string) to int
    try:  
        temp = float(temp)  
    except ValueError:  
        # temp was not a number, so silently
        # ignore/discard this line
        continue
    # this IF-switch only works because Hadoop sorts map output
    # by key (here: filename) before it is passed to the reducer
    if current_name == filename:  
        current_temp = max(temp, current_temp) 
    else:  
        if current_name:  
            # write result to STDOUT
            print('%s\t%s' % (filename, current_temp))
        current_temp = 0  
        current_name = filename  
# do not forget to output the last filename if needed!
if current_name == filename:  
    print('%s\t%s' % (current_name, current_temp))