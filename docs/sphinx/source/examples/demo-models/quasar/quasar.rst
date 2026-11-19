Demo: Quasar, M20
########################################################

The collaboration has published a series of papers using parameterised, biconical disc wind models. The initial model focus mostly on broad absorption quasars (`Higginbottom et al 2013 <https://ui.adsabs.harvard.edu/abs/2013MNRAS.436.1390H/abstract>`_), since the emission line were too weak in that case too match observed BLR properties. In `Matthews et al 2016 <https://ui.adsabs.harvard.edu/abs/2016MNRAS.458..293M/abstract>`_, we included a treatment of clumping and found some. Finally, in `Matthews et al 2020 <https://ui.adsabs.harvard.edu/abs/2020MNRAS.492.5540M/abstract>`_ (hereafter M20) we used a similar model to the previous clumpy wind model, but explored some of the behaviour in the ionizing flux density plane, and also used an isotropic illuminating SED.

This particular document focuses on Model A from M20. As with most of the demo models discussed here, the model makes use of the :doc:`Shlosman \& Vitello (1993) wind prescription<../../../wind_models/sv>`.

The wind is equatorial, and illuminated by an isotropic SED. 

.. note::

    more description needed

Important Parameters
============================
Central Source Parameters:

.. math::

	M_{\rm BH}				&= 	10^9 M_\odot			 	\\
	\dot{M}_{\rm acc}		&= 	5 M_\odot~{\rm yr}^{-1}	  	\\
	L_{2-10~{\rm keV}}		&=	10^{43}~{\rm erg~s}^{-1}			

Wind parameters:

.. math::

	\dot{M}_{\rm wind}  &=  5 M_\odot~{\rm yr}^{-1} \\
	\theta_{\rm min}	&= 	70^\circ \\	 
	\theta_{\rm max}	&= 	85^\circ \\	 
	r_{\rm min}     	&=	300 r_g \\		
	r_{\rm max}   		&= 	600 r_g \\	
	R_v 				&=  10^{19}~{\rm cm} \\
	f_V           		&=	0.01 

Illuminating SED 
============================

Runtime 
============================
3h25min on 64 cores (218 core hours).

Outputs 
============================

References
============================
