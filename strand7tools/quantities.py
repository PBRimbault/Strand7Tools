"""
The force quantities of a member force diagram, in the MEMBER axes, tension positive:

  N   axial force
  Vy  shear along y   (Strand7 SF1, along principal axis 1)
  Vz  shear along z   (Strand7 SF2, along principal axis 2: the web of an I-section)
  T   torque
  My  moment about y, MAJOR-axis bending   (Strand7 BM2)
  Mz  moment about z, minor-axis bending   (Strand7 BM1)

Diagrams are in kN and kNm; the EC3 modules work in N and N*m.
"""

QUANTITIES = ["N", "Vy", "Vz", "T", "My", "Mz"]
