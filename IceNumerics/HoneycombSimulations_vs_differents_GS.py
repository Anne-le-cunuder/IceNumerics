
# coding: utf-8

# # Honeycomb Colloidal Ice Simulations
# 
# This script makes simulations of Honeycomb Colloidal Ice, using LAMMPS as a Brownian Dynamics Engine. The setup of the input files for LAMMPS is done through a series of Python programs. The ouptut of the simulations can be either kept in the same file, or grouped together in another directory. 
# 
# This script uses several libraries that need to be imported. Some are "homemade" and some are external libraries. Make sure the external libraries are installed (most of them can be installed by pip)
# 
# The external libraries are:
# * numpy <- Numerical Arrays for python
# * matplotlib <- drawing plots for python
# 
# The homemade libraries are 
# * Spins <- Create Spin Ice Geometries.
# * ColloidalIce <- Turn Spins objects to Colloidal Ice by adding physical properties. 
# * LAMMPSInterface <- This takes Colloidal Ice objects and writes inputs scripts that LAMMPS can run.
# 

# In[ ]:


from Spins import *
from ColloidalIce import ColloidalIce
from LAMMPSInterface import *
import subprocess # Subprocess is a default library which allows us to call a command line program from within a python script
import shutil # shutil allows us to move files around. This is usefull to organize the resulting input and output files.
import matplotlib
#get_ipython().magic('matplotlib inline')


# Check that the imports work ok. If they don't, it means that you have to install something else.
# 
# Now we will use these tools to build a simulation setup. The command:

# In[ ]:


S = HoneycombSpinIce(10,10,Lattice = 30e3, Periodic = False, Ordering = "GroundState3");
fig1 = plt.figure(1);
ax1 = plt.axes();
for s in S:
    S[s].display(ax1)
plt.axis("equal");
plt.axis("tight");


# In[ ]:





# creates a Honeycomb geometry with 6x6 unit cells. 
# * *Lattice = 30e3* - is the lattice constant in nanometers (that is, the distance between two unit cells)
# * *Ordering = "Random"* - makes it a random initial configuration. You can also write *"Bias"* to get a bias initial state. 

# In[ ]:


Params = SimulationParameters(
    Runs=1,
    Thermo=1e2,
    Timestep= 0.5e-2,
    Framerate = 1,Time = 10)


# This command sets up the parameters used by the simulation. 
# 
# * *Runs=50* - This is the number of experiments performed in parallel. 
# * *Thermo=1e2* - This is how often (in units of timestep) LAMMPS prints output. It is better if *thermo* is not too small as to slow down the simulation, but enough that you can tell that it hasn't crashed (although LAMMPS rarelly crashes in this context)
# * *Timestep= 1e-2* - This is the size of the timestep in seconds. 
# * *Framerate = 10* - This is the frequency with which LAMMPS writes output (in seconds) a higher frequency will print more detailed output, but that output will be slower to read. Remember the output has the result of all the parallel simulations, so it can easily become very heavy.
# * *Time = 60* - This is the simulation duration (in seconds). 
#     

# In[ ]:


TargetDir = '.\\Output\\'


# *Target Dir* is a variable that determines where the results from the simulations are stored. 

# In[ ]:


C = ColloidalIce(S,FieldZ=[5,0],Stiffness_Spread = 0.2, height=200)
fig1 = plt.figure(1);
ax1 = plt.axes();
for c in C:
    C[c].display(ax1)
plt.axis("equal");
plt.axis("tight");


# The FieldZ parameter gives has two values. The first one is the strength of the maximum field. The second parameter is the time it takes to reach that field. 

# In[ ]:


L = LAMMPSScript(C,Params)
subprocess.call("lmp_mingw64.exe -in "+L.filename+".CI")
shutil.move(L.filename+".CI",TargetDir+"Data and Run Files\\"+L.filename+".CI")
shutil.move(L.filename+".data",TargetDir+"Data and Run Files\\"+L.filename+".data")
shutil.move(L.filename+".lammpstrj",TargetDir+L.filename+".lammpstrj")


# We will now attempt to lazy read the final state of the output file, and to visualize it first as a colloidal system, then as a spin system. 

# In[ ]:


Time = 600;
Exp = 1;
OutputFilename = TargetDir+L.filename+".lammpstrj"


# There are three relevant commands to read a data file. 
# * readline() reads the line that is after the line that was read last.
# * seek(byte) goes to the position *byte* in the file and then reads the rest of the line
# * tell() gives the current position on the file. 

# First we open the file

# In[ ]:


DataFile = open(OutputFilename,'r')


# Then we make a time dictionary. A time dictionary is a structure that stores the position in the file where I can find the timestep t (in units of timesteps). That is:
# * TimeDict[Time] <- gives the byte in the file where the data for the timestep Time starts
# 
# This is done by scanning the entire file and storing the value of "tell" when the keyword "TIMESTEP" is found

# In[ ]:



TimeDict={}
Line = DataFile.readline()
while Line:
    if 'TIMESTEP' in Line:
        Line = DataFile.readline()
        TimeDict[int(Line)] = DataFile.tell()
    else:
        Line = DataFile.readline()


# Now we can go to the byte where the requested timestep is located:

# In[ ]:


DataFile.seek(TimeDict[Time])


# Next we get rid of the information about number of atoms, and box bounds which is found at the start of every timestep:

# In[ ]:


Line = DataFile.readline()
while 'ITEM: ATOMS' not in Line:
    Line = DataFile.readline()


# Now we read line by line in the file. 
# 
# The standard .lammpstrj file format is "id type x y z". Particles in the same experiment have the same "type" property. To only store the requested experiment, we first convert the second element in the line (the function split(delim) splits a string into arrays of strings separated by the delimiter delim) to a float, and if it is equal to the requested experiment we convert all the line and store it in the list Row. 
# 
# In the end, we convert the list of lists Row to a numpy array

# In[ ]:


Row = []
while Line!="" and 'TIMESTEP' not in Line:
    Line = DataFile.readline()
    if Line!="" and 'TIMESTEP' not in Line:
        if float(Line[0:-2].split(' ')[1])==Exp:
            Row = Row+[[float(i) for i in Line[0:-2].split(' ')]]

Row = np.array(Row)


# In[ ]:


plt.plot(Row[:,2],Row[:,3],'+')


# In[ ]:


DataFile.close()


# In[ ]:




