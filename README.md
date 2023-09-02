export CPLEX_STUDIO_DIR="/opt/ibm/ILOG/CPLEX_Studio201/"
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/opt/ibm/ILOG/CPLEX_Studio201/opl/bin/x86-64_linux/
export JULIA_COPY_STACKS=1
julia --threads=auto --project=. data/campinas-random/1.sbrp --out solutions/complete_campinas-random1 --intersection-cuts --ip
