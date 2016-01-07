# NOTE: there is more simulation data contained in directory data

../../alchemical_analysis//alchemical_analysis.py \
  -a Gromacs \
  -m all \
  -d data/3-methylindole-11steps \
  -p dhdl \
  -q xvg \
  -o . \
  -r 5 \
  -u kcal \
  -s 0 \
  -c \
  -f 10 \
  -g \
  -w \
  > run.log
