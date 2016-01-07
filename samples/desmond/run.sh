../../alchemical_analysis//alchemical_analysis.py \
  -a Desmond \
  -m bar \
  -t 298.0 \
  -d . \
  -p 'gibbs.[0-9]*' \
  -q dE \
  -o . \
  -r 5 \
  -u kcal \
  -s 0 \
  -g \
  > run.log
