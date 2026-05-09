# Ko 2018 MEAM potential for Sn (NIST IPR id: 2018--Ko-W-S-Kim-D-H-Kwon-Y-J-Lee-M--Sn)

variable POTDIR string "<PROJECT_ROOT>/meam_potentials"

pair_style	meam
pair_coeff	* * ${POTDIR}/library.Sn.Ko2018.meam Sn ${POTDIR}/Sn.Ko2018.meam Sn

neighbor 1.0 bin
neigh_modify once no every 1 delay 0 check yes

min_style	     cg
min_modify	     dmax ${dmax} line quadratic

thermo		1
thermo_style custom step temp pe press pxx pyy pzz pxy pxz pyz lx ly lz vol
thermo_modify norm no
