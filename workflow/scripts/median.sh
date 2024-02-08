#  GNU nano 4.8                                                                                    median.sh                                                                                    Modificato  
sort -n $1 | awk '{
  data[NR] = $1
}
END {
  n = NR
  if (n % 2 == 0) {
    median = (data[n / 2] + data[n / 2 + 1]) / 2
  } else {
    median = data[(n + 1) / 2]
  }
  printf "Median: %.2f\n", median
}'



