SRoll method
=====

.. _methods:


Step[1]:  Reorganize data - Create Healpix Rings
------------------------
From Time Ordered Information, the algorithm create pixelized data called healpix rings (HPR) using HEALPIX library.


HEALPIX(Hierarchical Equal Area isoLatitude Pixelization)
**********

In the HEALPix standard the pixelization of the sphere begins by dividing the sphere into 12 base pixels of equal area, as shown in figure below.There are four base pixels around each pole and four around the equator. 
Each of these base pixels can be subdivided into four child pixels, and subsequently, each child pixel can be subdivided into four pixels of a lower hierarchy level. 
This process can continue until the desired resolution is achieved.

.. image:: /images/healpix.jpg
  :height: 200
  :align: center
 

|



See documentation for `Healpy <https://healpy.readthedocs.io/en/latest/>`_


Step[2] : Calcul differences in the same HEALPIX pixel
------------------------------------------------------

SRoll fits systematic effects, noise in 1/f and calibration using differences between mesures in the same healpix pixel. 
The Sroll method is based on the fact that there is multiples observations of the signal. During the course of the mission there will be several observations for a single healpix pixel, the aim is to minimise the differences between theses observations.

Below is an example of a SRoll map showing only the rings passing through the carina area (in red on the map). 
There are several rings passing through the same healpix pixel, so several observations for the same pixel (here schematized by Px_scan1,Px_scan2).

.. image:: /images/cross_rings_carina.png
  :height: 300
  :width: 600
  :align: center
  

The algorithm is based on the redundancy of the observations of the signal in pixel made by the same detector at different times, or by different detectors in the same pointing period, to determine the response of each bolometer. 
It use conjugate gradient method to minimize differences between pixels (example below scan1_px,scan2_px) and some of the systematics effects are also remove by using pre-calulated templates.

.. image:: /images/carina_pix_calc.png
  :height: 300
  :width: 500
  :align: center
  

Step[3] :  Clean data and create maps
--------------------------------------
Clean the data using fitted parameters and projects the signal to create cleaned
maps. A map Sp is created from clean data using projection matrix as follow :

.. math::

    S_{p}=\big{(}(A^T_{p}A_{p})^{-1} A^T_{p}\big{)}M_{p} \nonumber

where:  
 :math:`M_{p}` is the measured signal in a pixel;
 
 :math:`(A^T_{p}A_{p})^{-1} A^T_{p}\big{)}` is the projection matrix;
 
 
 The desired projection is define in parameters and define the output maps of SRoll. Hereafter an example of Sroll output project in Intensity and polarized (I,Q,U) :
 

.. image:: /images/IQUmaps.png
  :height: 300
  :width: 600
  :align: center
  
