#!/bin/zsh
cat << EOF
--- $0 --- Building LAMMPS
   (build) / clean
EOF

function DIE { echo "build error: $@">&2; exit -1; }

case `hostname` in
   R50) TARGET=gentoo ; ;;
   R52) TARGET=gentoo ; ;;
   server) TARGET=server; ;;
   lcmaster.rhrk.uni-kl.de) TARGET=rhrk ; ;;
   lcmaster1.rhrk.uni-kl.de) TARGET=rhrk; ;;
   earth) TARGET=servers; ;;
   wap) TARGET=wap; ;;
   login3) TARGET=hercules; ;;
   *) DIE "missing target arch" ; ;;
esac

function chk_src_dir {
   case `basename \`pwd\`` in 
      src) echo; ;;
      *) cd src ; ;;
   esac
}

function cleanup {
   chk_src_dir
   rm lmp_*
   make clean
   cd STUBS
   make clean 
   cd ..
}

function make_stubs {
   cd STUBS
   make clean
   make 
   cd ..
}

function make_packages {
   make yes-opt
}

function tag_version {
   TMPFILE=/tmp/tmp.lammps.`cat /dev/urandom | tr -cd a-z | head -c 5`
   cat universe.cpp | awk -v tag=`date +%d%m%Y-ziegen@rhrk.uni-kl.de` '/version = /{print "version = \""tag"\";"}!/version/{print}' > $TMPFILE
   mv $TMPFILE universe.cpp
}

function build {
   chk_src_dir
   make_stubs
   make_packages
   tag_version
   case $TARGET in
     gentoo)
        make -j 4 gentoo 2>&1| grep -v deprecated | grep -v "In member functio"
     ;;
     serial)
        make -j 4 serial
     ;;
     server)
        make -j 4 serial
        make -j 4 server_parallel
     ;;
     rhrk)
        make -j 4 rhrk_serial
        make -j 4 rhrk_parallel
     ;; 
     wap)
        make -j 4 wap_parallel
        make -j 4 serial
     ;;
     hercules)
        make -j 4 serial
        make -j 4 hercules_parallel
        make -j 4 hercules_mvapich
     ;;
   esac
}

case $1 in
   clean) cleanup; ;;
   *) build; ;;
esac

