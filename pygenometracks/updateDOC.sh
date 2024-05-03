# This bash script can be used to regenerate the tables used
# in possible parameters when parameters changed in tracksclasses.
# As well as updating the usage part in the README.md
if [ -e docs/content/tracks/auto ]; then
    rm -r docs/content/tracks/auto
fi
mkdir -p docs/content/tracks/auto
python pygenometracks/getAllDefaultsAndPossible.py 
# This generates 2 tables that will be included in the readthedocs documentation
# As well as a file in rst format for each track
echo "All available tracks and types
==============================

.. toctree::
   :maxdepth: 1
" > docs/content/all_tracks.rst
for f in docs/content/tracks/auto/*_deduced_from_code.txt; do
    track_name=`basename $f _deduced_from_code.txt`
    echo "   tracks/${track_name}" >> docs/content/all_tracks.rst
    if [ ! -e docs/content/tracks/${track_name}.rst ]; then
        echo "$track_name
==========

Description
-----------

Parameters
----------

.. include:: auto/${track_name}_deduced_from_code.txt

Output of \`\`make_tracks_file\`\`:
-------------------------------

.. literalinclude:: auto/${track_name}_options_text.txt
    :language: INI

">docs/content/tracks/${track_name}.rst
    fi
done

# Then for the usage:
awk 'NR==1,/<!--- Start of possible arguments of pgt -->/' README.md > newREADME.md
echo "\`\`\` text" >> newREADME.md
pyGenomeTracks -h | awk '/optional arguments:/{toprint = 1}toprint == 1 {print}' >> newREADME.md
echo "\`\`\`" >> newREADME.md
awk '/<!--- End of possible arguments of pgt -->/{toprint = 1}toprint == 1{print}' README.md >> newREADME.md
mv newREADME.md README.md
