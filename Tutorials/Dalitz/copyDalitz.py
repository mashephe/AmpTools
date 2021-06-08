#!/usr/bin/env python
import sys
import os

if len(sys.argv) != 2:
  print "Usage:"
  print "  copyDalitz.py <newDirectory/newName>"
  print "     e.g. copyDalitz.py $HOME/PWA/Psi3Body"
  print "            creates a new example package called Psi3Body"
  print "            in the existing directory $HOME/PWA/"
  exit()

outputInfo = sys.argv[1];
for i in range(5):
  outputInfo = outputInfo.replace("//","/")
outputDirBase = outputInfo.rpartition("/")[0];
outputNameBase = outputInfo.rpartition("/")[2];
print "----------------------------"
print "-- creating a new package --"
print "----------------------------"
print "  directory: ", outputDirBase;
print "   new name: ", outputNameBase;
print ""
yn = str(raw_input("continue? (y/n) ")).upper()
if yn != "Y":
  exit()
print ""

for root, dirs, files in os.walk("."):
  if root == ".":
    root = ""
  if root[:2] == "./":
    root = root.replace("./","",1)
  inputDirBase = root
  root = root.replace("Dalitz",outputNameBase);
  root = root.replace("DALITZ",outputNameBase.upper());
  root = root.replace("dalitz",outputNameBase.lower());
  outputDirectory = outputDirBase+"/"+outputNameBase+"/"+root
  if root == "":
    outputDirectory = outputDirBase+"/"+outputNameBase
  if (root == "doc" or root == "lib"):
    continue
  print "make new directory: ", outputDirectory
  os.mkdir(outputDirectory)
  for name in files:
    if name[-3:] == ".py":
      continue
    inputName = inputDirBase+"/"+name
    if inputDirBase == "":
      inputName = name
    name = name.replace("Dalitz",outputNameBase);
    name = name.replace("DALITZ",outputNameBase.upper());
    name = name.replace("dalitz",outputNameBase.lower());
    outputName = outputDirectory+"/"+name
    print "make new file: ", outputName
    inputFile = open(inputName,"r")
    outputFile = open(outputName,"w")
    for line in inputFile:
      line = line.replace("Dalitz",outputNameBase);
      line = line.replace("DALITZ",outputNameBase.upper());
      line = line.replace("dalitz",outputNameBase.lower());
      outputFile.write(line)
    inputFile.close()
    outputFile.close()

