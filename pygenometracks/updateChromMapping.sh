if [ ! -e ChromosomeMappings ]; then
    git clone https://github.com/dpryan79/ChromosomeMappings.git
else:
    cd ChromosomeMappings
    git pull
    cd ..
fi
python pygenometracks/updateChromMapping.py
