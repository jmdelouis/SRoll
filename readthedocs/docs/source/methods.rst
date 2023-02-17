SRoll methods
=====

.. _methods:


Step[1]:  Reorganize data - Create Healpix Rings
------------------------
From fits data, Time Ordered Information, the algorithm create pixelized data called healpix rings (HPR) using HEALPIX library.


HEALPIX(Hierarchical Equal Area isoLatitude Pixelization)
**********

In the HEALPix standard the pixelization of the sphere begins by dividing the sphere into 12 base pixels of equal area, as shown in Figure 1(a).There are four base pixels around each pole and four around the equator. 
Each of these base pixels can be subdivided into four child pixels, and subsequently, each child pixel can be subdivided into four pixels of a lower hierarchy level. 
This process can continue until the desired resolution is achieved.

.. image:: /images/healpix.jpg
  :height: 150
  :align: center
 

|



See documentation for `Healpy <https://healpy.readthedocs.io/en/latest/>`_


Step[2] : Calcul differences in the same HEALPIX pixel
------------------------------------------------------

SRoll fits systematic effects, noise in 1/f and calibration using differences between mesures in the same healpix pixel.( see figure below). Some of the systematics effects are also remove by using pre-calulated templates.

.. image:: /images/scns.png
  :height: 300
  :width: 500
  :align: center
  
|

During the mission, the satellite observed the same area of the sky several times with a time shift of one rotation. In the case of the CMB observation the signal observed between a time t and t1(+1 rotation) does not change or changes very little, so the hypothesis is 
that the differences between theses observations are instrumental  or foreground effects . We then extract these offsets by calculating the difference between observation at t and t1.

Step[3] :  Clean data and create maps
--------------------------------------
Clean the data using fitted parameters and projects the signal to create cleaned
maps. A map Sp is created from clean data using projection matrix as follow :

.. math::

    S_{p}=\big{(}(A^T_{p}A_{p})^{-1} A^T_{p}\big{)}M_{p} \nonumber

where:  
 :math:`M_{p}` is the measured signal in a pixel;
 
 :math:`(A^T_{p}A_{p})^{-1} A^T_{p}\big{)}` is the projection matrix;
 
 
SRoll generic data model
--------------------------------------

.. math::

g_{d}M_{d,r,p,t} = A_{d,r,p,t} \cdot S_{p,t} + \sum_{h=1}^{n_{\mathrm{sys}}}\gamma_{d,h} T_{d,r,p,t,h} + 
\sum_{f=1}^{n_{\mathrm{comp}}} L_{d,f} C_{p,t,f} + O_{d,r} + g_{d}N_{d,r,p,t},

where:

 :math:` $M_{d,r,p,t}$` is the measured detector total signal of a detector chain $d$, during the stationary period $r$, inside the pixel $p$ at the date $t$;
 :math:`$g_{d}$` is the calibrated response of a detector chain $d$;
\item $A_{d,r,p}$  is the pointing vector giving the observed pixel for a given detector in a given ring (e.g. $\left[1,\rho_d\cos(2\phi_{d,r,p}),\rho_d\sin(2\phi_{d,r,p})\right]$ for Planck polarized map; $\rho_b$ is the ground-measured polarization efficiency;$\phi_{d,r,p}$ is the ground-measured detector polarization angle for detector $d$ with respect to the north-south axis;);

\item $S_{p,t}$ is the sky signal in pixel $p$ after subtraction of the orbital dipole assumed to be known with an amplitude invariant in time (e.g  $\left[I_p, Q_p, U_p\right]$ for Planck polarized map where $I_p$, $Q_p$, and $U_p$ represent the common sky maps seen by all detectors).;

\item $T_{d,r,p,t,h}$ is the spatiau temporal template of an instrumental systematic effect. This template can be gicen to sroll or synthetised through generative neural network using an internal loop.

\item $\gamma_{b,h}$ is the amplitude of this systematic effect;

\item $C_{p,t,f}$ are the foreground component spatiau temporel templates which can be given in input or synthetise within an internal loop using generative neural network;
\item $L_{d,f}$ is the detector efficient for the $f$ foreground;
\item $O_{d,r}$ is the offset per pointing period $r$ used to model the $1/f$ noise, and we set $\sum_{b=1}^{n_{\mathrm{bolo}}} \sum_{r=1}^{n_{\mathrm{ring}}}O_{b,i}=0$, since the SRoll algorithme is based on differences and can not extract information on the monopole;
\item $N_{d,r,t,p}$ is the pixel white noise (electronic and/or photon noises), with variance ${\sigma_{d,r,t,p}}^2$;
\end{itemize}
It is important to understand that the purpose of the SRoll algorithm is not only to measure $S_{p,t}$, but also to extract the best possible models from foregrounds $C_{p,t,f}$ and instrumental systematic effects $T_{d,r,p,t,h}$. Also, the signal of scientific interest can be part of the $S_{p,t}$, or the foregrounds $C_{p,t,f}$. This notion of the signal of interest depends on the data on which Sroll is applied.
