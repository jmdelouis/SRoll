Documentation - SRoll4
===================================

**SRoll** is algorithm based on data inversion and
scattering transform to reduce systematics effects and foregrounds on satellite data. You can find the evolutions of the SRoll algorithm and the distribution of its products here : `sroll2`_

.. _sroll2: http://sroll20.ias.u-psud.fr


Context 
********

Data from space observations have instrumental effects or foregrounds that degrade the signal of interest recovery. The SRoll software developed initialy for
cosmology (`PLANCK`_) simultaneously inverts this signal and the instrumental effects including noises. Therefore, SRoll uses the all available data set
to minimise the variance of the measured calibration parameters. This problem is common to astrophysical and oceanographic data as with `CFOSAT`_, but the
geophysical signal evolves over time on the scale of the measurement, unlike astrophysical processes. Therefore, we present a new version of SRoll based on
Scattering Transforms (ST) to model the signal dynamics. STs have been used successfully on interstellar medium thanks to the intermittent structure of its
turbulent processes. By using these statistical constraints measured on the data itself there is no learning phase as with usual neural networks.


.. _CFOSAT: https://cfosat.cnes.fr/fr/CFOSAT/Fr/index.htm

.. _PLANCK: https://www.esa.int/Enabling_Support/Operations/Planck

Citations
*********

[1] Vibert L. Delouis Jean Marc Puget J.-L. Improved large-scale interstellar dust foreground model and CMB solar dipole measurement. Ed. by EDP
Sciences. Article. FRANCE, 2021. doi: https://doi.org/10.1051/0004-6361/202140616. url: https://archimer.ifremer.fr/doc/00702/81365/.

[2] K. M. Gorski et al. “HEALPix: A Framework for High-Resolution Discretization and Fast Analysis of Data Distributed on the Sphere”. In: The
Astrophysical Journal 622.2 (Apr. 2005), pp. 759–771. doi: 10.1086/427976. url: https://doi.org/10.1086%2F427976.

[3] M. Lopez-Radcenco, J.-M. Delouis, and L. Vibert. “SRoll3: A neural network approach to reduce large-scale systematic effects in the iPlanck/i
High-Frequency Instrument maps”. In: Astronomy &amp Astrophysics 651(July 2021), A65. doi: 10.1051/0004-6361/202040152. url: https://doi.org/10.1051%2F0004-6361%2F202040152.

[4] Pagano, L. et al. “Reionization optical depth determination from Planck HFI data with ten percent accuracy”. In: A&A 635 (2020), A99. doi:10.1051/0004-6361/201936630. url: https://doi.org/10.1051/0004-6361/201936630.

[5] Delouis, J. M., Allys, E., Gauvrit, E., & Boulanger, F. (2022). Non-Gaussian modelling and statistical denoising of Planck dust polarisation full-sky maps using scattering transforms. Astronomy & Astrophysics, 668, A122.DOI: https://doi.org/10.1051/0004-6361/202244566 

.. note::

   This project is funded by a CNES (National Center for Space Studies) R&T.

Contents
--------

.. toctree::

   installation
   get_started
   list_parameters
   methods
