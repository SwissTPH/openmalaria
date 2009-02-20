#!/bin/sh
rm test/sandbox/* 2>/dev/null
cp model/openMalaria test/sandbox
cd test/sandbox && cp ../original/* . 2>/dev/null
cp scenario2.xml scenario.xml
echo "Running scenario2.xml to create init_data.xml"
./openMalaria  > /dev/null
wait
if [ -e checkpoint ] 
  then
  rm checkpoint*
fi

#set the checkpointing period to 10 secs
sed -e 's/<checkpoint_period>.*<\/checkpoint_period>/<checkpoint_period>10<\/checkpoint_period>/' init_data.xml > init_data_checkpoint.xml
mv init_data_checkpoint.xml init_data.xml
cp checkpointTest.xml scenario.xml
echo "Running checkpointTest.xml to create checkpoint files"
./openMalaria  > /dev/null 
wait
#test to make sure there is a checkpoint file
if [ -e checkpoint ]
  then 
  rm output.txt
  echo "Running checkpointTest.xml from checkpoint"
  ./openMalaria  > /dev/null
  wait
  test -e output.txt && ../original/compareOutputsFloat.py originalCheckpointTest.txt output.txt 1
else
echo "No checkpoint file found, testing of checkpointing not possible"
  exit
fi

