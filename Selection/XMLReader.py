import xml.etree.ElementTree as ET

ifilename='selection2016.xml'

#tree = ET.parse(ifilename)
tree = ET.parse('selection2016.xml')
#tree = ET.parse('test.xml')
root = tree.getroot()

#root = ET.fromstring(country_data_as_string)

for child in root:
 print child.tag, child.attrib

import SelectionMaker as sel
s = sel.Selection()
for var in root.iter("Variable"):
   print var.attrib
   s.AddVariable(var.attrib)


for bins in root.iter("Bin"):
   print bins.attrib
   print type(bins.attrib)
   s.AddBin(bins.attrib)


s.CreateSelFunctions()

#for neighbor in root.iter('neighbor'):
#	print neighbor.attrib

#for v in root.iter("var"):
#   print v

for v in root.findall('var'):
	print "toto"
	print v.get('name')

for country in root.findall('country'):
    rank = country.find('rank').text
    name = country.get('name')
    print name, rank




