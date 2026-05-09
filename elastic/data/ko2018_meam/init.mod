# beta-Sn (white tin, A5 prototype, I4_1/amd) for Ko 2018 MEAM
# 4-atom conventional tetragonal cell

variable up equal 1.0e-4
variable atomjiggle equal 1.0e-5

units		metal
variable cfac equal 1.0e-4
variable cunits string GPa

variable etol equal 0.0
variable ftol equal 1.0e-10
variable maxiter equal 200
variable maxeval equal 2000
variable dmax equal 1.0e-2

# Initial guess from NIST IPR (Ko 2018 predictions)
variable a equal 5.8589
variable c equal 3.2059

boundary	p p p

lattice custom 1.0 &
  a1 ${a} 0.0 0.0 &
  a2 0.0 ${a} 0.0 &
  a3 0.0 0.0 ${c} &
  basis 0.0 0.0 0.0 &
  basis 0.5 0.5 0.5 &
  basis 0.0 0.5 0.25 &
  basis 0.5 0.0 0.75

region		box prism 0 1 0 1 0 1 0.0 0.0 0.0
create_box	1 box
create_atoms	1 box

mass 1 118.71
