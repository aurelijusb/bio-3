Bioinfomratics
==============

University project for Bioinformatics.

Compiling/running
-----------------

    python HpvAnalyser.py
    
Libraries and external tools
----------------------------

 * BioPython http://biopython.org/wiki/Main_Page
 * CD-HIT http://weizhong-lab.ucsd.edu/cd-hit/
 * MAfft http://mafft.cbrc.jp/alignment/software/
 * matplotlib http://matplotlib.org/
 * gtk+ http://www.gtk.org/
 * cairo http://www.cairographics.org/

Not finished
------------

Finding similar harmful probes and disimilar safe probes.
No elements with needed length. Probably errors in earlyer algorythm.
GTK and cairo used only for debuging purposes.

Current version is not stable. Use earlier revisions for more stable
(but not finished) application versions.

Functionality
-------------

 * Generates query to get all Human papillomavirus types
 * Retrieve data from nucleotide database
 * Use CD-HIT to remove duplicates
 * Use Mafft to align sequences
 * Saves results and caching files (per each HPV type, unaligned and aligned)
 
Not fully finised:

 * Safe conservation by each base from aligned
 * Visualise aligned sequences using gtk+ and cairo
 * Use brute-force to calculate simlilar harmful probes (with error allowed)

References
----------

 * http://biopython.org/DIST/docs/tutorial/Tutorial.html
 * http://bips.u-strasbg.fr/fr/Tutorials/Comparison/Blast/blastall.html
 * http://ghr.nlm.nih.gov/glossary=probe
 * http://www.ncbi.nlm.nih.gov/nuccore/
