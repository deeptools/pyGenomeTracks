# This bash script can be used to regenerate the README when parameters changed in tracksclasses and/or ini files for alpha and narrowpeak:

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

# Include the ini file of alpha:
awk 'NR==1,/<!-- Master alpha ini start-->/' README.md > newREADME.md
echo "\`\`\`INI" >> newREADME.md
cat pygenometracks/tests/test_data/alpha.ini >> newREADME.md
echo "\`\`\`" >> newREADME.md
awk '/<!-- Master alpha ini end-->/{toprint = 1}toprint == 1{print}' README.md >> newREADME.md
mv newREADME.md README.md


# Include the ini file of narrowPeak:
awk 'NR==1,/<!-- Master narrow ini start-->/' README.md > newREADME.md
echo "\`\`\`INI" >> newREADME.md
cat pygenometracks/tests/test_data/narrow_peak.ini >> newREADME.md
echo "\`\`\`" >> newREADME.md
awk '/<!-- Master narrow ini end-->/{toprint = 1}toprint == 1{print}' README.md >> newREADME.md
mv newREADME.md README.md
