Line_transfer
=============
The way in which line transfer and scattering is dealt with
in the code. Governs whether we adopt any approximations
for radiative transfer, whether to use the indivisible packet
and macro-atom machinery, and whether to use isotropic or
anisotropic scattering.

The classic mode is recommended for non-macro atom runs,
while the macro mode is recommended for macro-atom runs.

The older names for these modes (escape_prob, thermal_trapping,
macro_atoms_escape_prob and macro_atoms_thermal_trapping) are still accepted as
synonyms for classic_iso, classic, macro_iso and macro respectively, but
are not shown in the prompt and are replaced by the new names in the .out.pf file.

Type
  Enumerator

Values
  pure_abs
    *Pure absorption*
    
    The pure absortion approximation.

  pure_scat
    *Pure scattering*
    
    The pure scattering approximation.

  sing_scat
    *Single scattering*
    
    The single scattering approximation.

  classic_iso
    *Escape probability + isotropic scattering* (synonym: escape_prob)
    
    Resonance scattering and electron scattering is dealt with isotropically.
    free-free, compton and bound-free opacities attenuate the weight of the photon
    wind emission produces additional photons, which have their directions chosen isotropically.
    The amount of radiation produced is attenuated by the escape probability.

  classic
    *Escape probability + anisotropic scattering* (synonym: thermal_trapping)
    
    As classic_iso, but we use the 'thermal trapping method' to choose an
    anistropic direction when an r-packet deactivation
    or scatter occurs.

  macro_iso
    *Macro-atoms + isotropic scattering* (synonym: macro_atoms_escape_prob)
    
    use macro-atom line transfer.
    Packets are indivisible and thus all opacities are dealt with by activate a macro-atom, scattering,
    or creating a k-packet.
    the new direction following electron scattering or deactivation of
    a macro atom is chosen isotropically.

  macro
    *Macro-atoms + anisotropic scattering* (synonym: macro_atoms_thermal_trapping)
    
    as macro_iso, but we use the 'thermal trapping method' to choose an anistropic direction
    when an r-packet deactivation or scatter occurs.


File
  `setup_line_transfer.c <https://github.com/sirocco-rt/sirocco/blob/master/source/setup_line_transfer.c>`_


Child(ren)
  * :ref:`Reverb.matom_lines`

  * :ref:`Wind_heating.kpacket_frac`

