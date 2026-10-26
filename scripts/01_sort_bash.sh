# Problem 1: run from IBS7048_2026 (the folder that contains raw/)
for f in raw/*.fastq
do
  condition=$(echo "$f" | cut -d '_' -f 2)
  r_direction=$(echo "$f" | cut -d '_' -f 4 | cut -d '.' -f 1)
  mkdir -p "sorted_bash/$condition/$r_direction"
  cp "$f" "sorted_bash/$condition/$r_direction/"
done
echo "raw: $(ls raw | wc -l) files"
for d in sorted_bash/*/*/
do
echo "$d: $(ls "$d" | wc -l) files"
done
