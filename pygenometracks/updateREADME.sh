# This bash script can be used to regenerate the tables used
# in possible parameters when parameters changed in tracksclasses.
# As well as updating the usage part in the README.md

python pygenometracks/getAllDefaultsAndPossible.py 

# This generates 2 tables that will be included in the readthedocs documentation

# Then for the usage:
awk 'NR==1,/<!--- Start of possible arguments of pgt -->/' README.md > newREADME.md
echo "\`\`\` text" >> newREADME.md
bin/pyGenomeTracks -h | awk '/optional arguments:/{toprint = 1}toprint == 1 {print}' >> newREADME.md
echo "\`\`\`" >> newREADME.md
awk '/<!--- End of possible arguments of pgt -->/{toprint = 1}toprint == 1{print}' README.md >> newREADME.md
mv newREADME.md README.md
