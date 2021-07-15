#!/bin/sh

PATH_TO_HADOOPSTREAM=/home/hadoop/hadoop-3.3.1/share/hadoop/tools/lib/hadoop-streaming-3.3.1.jar
PATH_TO_REDUCER=/home/hadoop/project3/reducer.py
PATH_TO_MAPPER=/home/hadoop/project3/mapper.py
PATH_OUTPUT=/tmp/output
echo "Running Hadoop Map Reduce..."

hadoop fs -rm -r $PATH_OUTPUT

hadoop jar $PATH_TO_HADOOPSTREAM \
-file $PATH_TO_MAPPER -mapper "mapper.py" \
-file $PATH_TO_REDUCER -reducer "reducer.py" \
-input /tmp/whether_data/* -output $PATH_OUTPUT