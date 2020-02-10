Usage
=====

.. contents:: 
    :local:

Starting usage
--------------

To run pyGenomeTracks a configuration file describing the tracks is required.
The easiest way to create this file is using the program ``make_tracks_file`` which creates a configuration file with defaults that can be easily changed.
``make_tracks_file`` uses the file ending to guess the file type. Then, a region can be plotted using ``pyGenomeTracks``. Both programs are described above:

make_tracks_file
----------------

.. argparse::
   :ref: pygenometracks.makeTracksFile.parse_arguments
   :prog: make_tracks_file
   :nodefault:

   
   
pyGenomeTracks
--------------

.. argparse::
   :ref: pygenometracks.plotTracks.parse_arguments
   :prog: pyGenomeTracks
   :nodefault:
