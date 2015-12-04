../../alchemical_analysis//alchemical_analysis.py \
  -a AMBER \
  -d data \
  -p '[01].*/ti00[2-9]' \
  -q out \
  -o . \
  -r 5 \
  -u kcal \
  -s 0 \
  -g \
  -w \
  > run.log
