#!/bin/bash
for i in *.R
do
  svn add $i
done
