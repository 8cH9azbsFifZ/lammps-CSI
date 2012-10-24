# Install/unInstall package classes in LAMMPS

if ($1 == 1) then

  cp style_colloid.h ..

  cp pair_colloid.cpp ..

  cp pair_colloid.h ..

else if ($1 == 0) then

  rm ../style_colloid.h
  touch ../style_colloid.h

  rm ../pair_colloid.cpp

  rm ../pair_colloid.h

endif
