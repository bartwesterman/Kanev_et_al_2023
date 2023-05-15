#!/bin/bash

INPUT_FILE=$1
java -jar ambit-tautomers-2.0.0-SNAPSHOT.jar -f $INPUT_FILE -o ${INPUT_FILE}_tautomer.csv -t best