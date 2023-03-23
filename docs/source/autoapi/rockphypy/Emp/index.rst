:py:mod:`rockphypy.Emp`
=======================

.. py:module:: rockphypy.Emp


Module Contents
---------------

Classes
~~~~~~~

.. autoapisummary::

   rockphypy.Emp.Empirical




.. py:class:: Empirical

   
   Empirical relations that widely applied
















   ..
       !! processed by numpydoc !!
   .. py:method:: esti_VS(vsh)

      
      Estimate, using the Greenberg-Castagna empirical relations, the shearwave velocity in a brine-saturated shaly sandstone with vp
      since we only assume two minearl phases: so L=2, X1= 1-vsh, X2= vsh

      :param Vp: compressional velocities, m/s
      :type Vp: float or array-like
      :param vsh: shale volume
      :type vsh: float or array-like

      :returns: * *float or array-like* -- vs (m/s): estimated shear wave velocities
                * *References* -- ----------
                * *- handbook of rock physics P516*















      ..
          !! processed by numpydoc !!

   .. py:method:: han(C)

      
      Han (1986) found empirical regressions relating ultrasonic (laboratory) velocities to porosity and clay content.effective pressure is 20Mpa

      :param phi: porosity
      :type phi: float or array-like
      :param C: clay volume fraction
      :type C: float or array-like

      :returns: *float or array-like* -- P and S wave velocities















      ..
          !! processed by numpydoc !!


