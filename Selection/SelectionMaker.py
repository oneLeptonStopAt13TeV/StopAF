
signs_ = ['>','>=','<',"<=",'==']

class Cut:
  
  def __init__(self):
     self.label = ""
     self.babytuple_name = ""
     self.cut_value = 0
     self.sign = ""


  def load(label, babytupleName, selection):
     self.label = label
     self.babytuple_name = babytupleName



  #def CreateSelFunction():


  #def AddSelection():



class Selection:


   def __init__(self):
	self.variables = []
	self.bins = []

   def AddVariable(self, variable_from_xml):
	self.variables.append(variable_from_xml)

   def AddBin(self, bin_from_xml):
	# change the selection with the following correspondance
	# and -> &&
	# !>  -> <
	bin = bin_from_xml
	bin['selection'] = bin['selection'].replace('and','&&', 100) # max nof occurence = 100
	bin['selection'] = bin['selection'].replace('!>','<', 100) # max nof occurence = 100
   	self.bins.append(bin)


   def CreateSelFunctions(self):
   	
	## replace MET by its babytuple name
	METbabytuple = "MET"
	for var in self.variables:
	  if var['name'] == "MET":
	    METbabytuple = "myEvent."+var['babyTuple']
	####################################

	## loop over the bins

	for bin in self.bins:
	   METbins = bin['METbins'].split(",")
	   
	   ## loop over the MET bins
	   for i  in range(len(METbins)-1):
	     #if i!=len(METbins)-1:
	       name=bin['name']+"_MET"+METbins[i]
	       selection = bin['selection']+" && "+METbabytuple+">="+METbins[i]
	       if METbins[i+1]!="inf" :
	         selection+=" && "+METbabytuple+"<"+METbins[i+1]  
		 name+="to"+METbins[i+1]
	       else: 
	         name+="toInf"
	       
	       #replace the variable label by their babytuple name
	       for var in self.variables:
	         selection = selection.replace(var['name'],"myEvent."+var['babyTuple'])
	      	 
	       dump="bool "+name+"() { "+selection+";}"
	       print dump

   #def AddSelection


