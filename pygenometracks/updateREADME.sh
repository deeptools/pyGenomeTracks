# This bash script can be used to regenerate the README when parameters changed in tracksclasses.

# Firts run the python script getAllDefaultsAndPossible.py
python pygenometracks/getAllDefaultsAndPossible.py 

# This generates 2 tables to include in the README
# First include the first one: docs/content/all_default_properties.txt
awk 'NR==1,/<!--- Start of default table -->/' README.md > newREADME.md
cat docs/content/all_default_properties.txt >> newREADME.md
awk '/<!--- End of default table -->/{toprint = 1}toprint == 1{print}' README.md >> newREADME.md
# Then include the second one: docs/content/all_possible_properties.txt
awk 'NR==1,/<!--- Start of possible table -->/' newREADME.md > README.md
cat docs/content/all_possible_properties.txt >> README.md
awk '/<!--- End of possible table -->/{toprint = 1}toprint == 1{print}' newREADME.md >> README.md
rm newREADME.md
