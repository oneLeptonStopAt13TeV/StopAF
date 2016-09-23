import xml.etree.ElementTree as ET
import os

###################################
# Goal is to create selection.h
###################################

###################################
# should become a parameter
###################################
ifilename='selection2016.xml'
template='template.h'
ofilename='test.h'

tree = ET.parse(ifilename)
#tree = ET.parse('selection2016.xml')
#tree = ET.parse('test.xml')
root = tree.getroot()


import SelectionMaker as sel
s = sel.Selection()

###################################
# Add Variable
###################################
for var in root.iter("Variable"):
   #print var.attrib
   s.AddVariable(var.attrib)

###################################
# Add Bins
###################################
for bins in root.iter("Bin"):
   #print bins.attrib
   s.AddBin(bins.attrib)

###################################
# Add Region
###################################
for regions in root.iter("Region"):
   s.AddRegion(regions.attrib)

###################################
# Create selection region
###################################
template='template.h'
ofilename='test.h'
ofilename2="test2.cc"
command = "cp "+template+" "+ofilename
os.system(command)

s.CreateSelFunctions(ofilename)
s.AddSelection(ofilename2)





